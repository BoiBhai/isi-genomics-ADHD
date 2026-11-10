# ============================================================
# Numeric-only MICE  ->  Maximal Elastic-Net with MI
# Outputs to: E:/gwas project/
# ============================================================
need <- c("readr","dplyr","stringr","tidyr","purrr","tibble",
          "mice","naniar","glmnet","Matrix","miceadds")
to_install <- setdiff(need, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(need, library, character.only = TRUE))

# ---------- CONFIG ----------
base_dir <- "E:/gwas project"
raw_csv  <- file.path(base_dir, "casecase.csv")
imputed_long <- file.path(base_dir, "casecase_numeric_imputed_long_exclGD.csv")

# Modeling config
Y_POST <- "pharma.post.ADHD Index"   # response
Y_PRE  <- "pharma.pre.ADHD Index"    # baseline covariate
DRUG   <- "pharma.Drug"

USE_MISSING_INDICATORS <- TRUE       # keep rows; treat missing genotypes as "UNK"
GENO_MIN_NONMISS_FRAC  <- 0.50       # include genotype markers with ≥50% observed (pre-UNK)
ALPHA_GRID <- seq(0.0, 1.0, by = 0.1)
N_FOLDS    <- 10
SEED       <- 20251105

DO_STABLE_POOLS  <- TRUE
STABILITY_THRESH <- 0.60             # selection frequency threshold for pooled OLS

# ---------- (A) Create numeric-only MICE long file if missing ----------
if (!file.exists(imputed_long)) {
  if (!file.exists(raw_csv)) stop("Raw CSV not found: ", raw_csv)
  message("Imputed long file not found. Creating numeric-only MICE from: ", raw_csv)
  
  df <- readr::read_csv(raw_csv, show_col_types = FALSE)
  
  # candidate numeric blocks to impute (ONLY these)
  num_blocks <- c("^expressions\\.", "^concentrations\\.", "^cprs-r\\.",
                  "^pharma\\.pre\\.", "^pharma\\.post\\.")
  numeric_cols <- unique(unlist(lapply(num_blocks, function(p) grep(p, names(df), value = TRUE))))
  
  # columns to exclude entirely from MICE (no impute, no predictor)
  genotype_cols <- grep("^genotype\\.", names(df), value = TRUE)  # includes genotype.Sex + markers
  id_col   <- intersect("ID", names(df))
  drug_col <- intersect("pharma.Drug", names(df))
  type_col <- intersect("Type", names(df))
  
  # coerce ONLY numeric candidates to numeric (observed values preserved)
  numify <- function(x) suppressWarnings(as.numeric(x))
  for (nm in numeric_cols) if (nm %in% names(df)) df[[nm]] <- numify(df[[nm]])
  
  dat_imp <- df
  meth <- make.method(dat_imp)
  pred <- make.predictorMatrix(dat_imp)
  
  # never impute non-numeric columns
  meth[setdiff(names(dat_imp), numeric_cols)] <- ""
  
  # numeric columns: PMM if at least 2 observed values; else skip
  for (nm in numeric_cols) {
    obs <- sum(!is.na(dat_imp[[nm]]))
    meth[nm] <- if (obs >= 2) "pmm" else ""
  }
  
  # exclude genotype.* + Drug + Type + ID from predictors/targets
  exclude <- unique(c(genotype_cols, drug_col, type_col, id_col))
  if (length(exclude)) {
    pred[, exclude] <- 0
    pred[exclude, ] <- 0
    meth[exclude] <- ""
  }
  
  # drop all-NA columns from prediction; no self-prediction
  all_na <- names(dat_imp)[colSums(!is.na(dat_imp)) == 0]
  if (length(all_na)) { pred[, all_na] <- 0; pred[all_na, ] <- 0; meth[all_na] <- "" }
  diag(pred) <- 0
  
  set.seed(SEED)
  imp <- mice(dat_imp, m = 20, maxit = 20, method = meth, predictorMatrix = pred, printFlag = TRUE)
  
  readr::write_csv(complete(imp, action = "long", include = TRUE), imputed_long)
  message("Saved imputed long file: ", imputed_long)
} else {
  message("Found imputed long file: ", normalizePath(imputed_long, winslash = "/"))
}

# ---------- (B) MAXIMAL ELASTIC-NET ----------
Dlong <- readr::read_csv(imputed_long, show_col_types = FALSE)
if (!(".imp" %in% names(Dlong))) stop("Stacked long imputations (.imp) not found in: ", imputed_long)
Dlong <- dplyr::filter(Dlong, .imp >= 1)
m <- max(Dlong$.imp, na.rm = TRUE)

# Feature families present
all_names <- names(Dlong)
geno_cols <- grep("^genotype\\.M\\d+", all_names, value = TRUE)  # M01..M22 markers
num_cols  <- unique(unlist(lapply(
  c("^expressions\\.", "^concentrations\\.", "^cprs-r\\.", "^pharma\\.pre\\."),
  function(p) grep(p, all_names, value = TRUE)
)))
num_cols <- setdiff(num_cols, c(Y_POST, Y_PRE))  # avoid duplicate baseline

# existence checks
if (!all(c(Y_POST, Y_PRE, DRUG) %in% all_names)) {
  stop("Missing required columns: ", paste(setdiff(c(Y_POST, Y_PRE, DRUG), all_names), collapse=", "))
}

# genotype inclusion filter (observed fraction BEFORE adding UNK)
pick <- function(v) intersect(v, all_names)
geno_cols <- pick(geno_cols)
num_cols  <- pick(num_cols)

obs_frac <- if (length(geno_cols)) sapply(Dlong[geno_cols], function(x) mean(!is.na(x))) else numeric(0)
geno_keep <- if (length(obs_frac)) names(obs_frac[obs_frac >= GENO_MIN_NONMISS_FRAC]) else character(0)
message("Using ", length(geno_keep), " genotype markers (≥ ", 100*GENO_MIN_NONMISS_FRAC, "% observed).")

# per-imputation design builder
prep_one_imp <- function(d) {
  d <- d %>% dplyr::filter(!is.na(.data[[Y_POST]]), !is.na(.data[[Y_PRE]]), !is.na(.data[[DRUG]]))
  d[[DRUG]] <- factor(d[[DRUG]])
  
  if (length(geno_keep)) {
    for (g in geno_keep) {
      if (USE_MISSING_INDICATORS) {
        d[[g]] <- factor(ifelse(is.na(d[[g]]), "UNK", as.character(d[[g]])))
      } else {
        d <- d[!is.na(d[[g]]), , drop = FALSE]
        d[[g]] <- factor(d[[g]])
      }
    }
  }
  
  rhs_terms <- c(
    paste0("`", Y_PRE, "`"),
    paste0("factor(`", DRUG, "`)"),
    if (length(geno_keep)) paste0("factor(`", geno_keep, "`)") else NULL,
    if (length(geno_keep)) paste0("factor(`", DRUG, "`):factor(`", geno_keep, "`)") else NULL,
    if (length(num_cols))  paste0("`", num_cols, "`") else NULL
  )
  rhs <- paste(rhs_terms, collapse = " + ")
  form_mm <- stats::as.formula(paste0("~ 0 + ", rhs))   # no intercept; glmnet adds it
  
  X <- model.matrix(form_mm, data = d)
  if (ncol(X) > 0) {
    keep_var <- which(Matrix::colSums((X - Matrix::colMeans(X))^2) > 0)
    if (length(keep_var)) X <- X[, keep_var, drop = FALSE]
  }
  list(y = d[[Y_POST]], X = X, data = d)
}

Dlist <- split(Dlong, Dlong$.imp)
designs <- lapply(Dlist, prep_one_imp)

# union of features across imputations; align each X to this set
all_feats <- Reduce(union, lapply(designs, function(z) colnames(z$X)))
align_X <- function(X, feats) {
  if (is.null(colnames(X))) {
    X <- matrix(0, nrow = nrow(X), ncol = length(feats)); colnames(X) <- feats; return(X)
  }
  miss <- setdiff(feats, colnames(X))
  if (length(miss)) {
    add <- matrix(0, nrow = nrow(X), ncol = length(miss)); colnames(add) <- miss
    X <- cbind(X, add)
  }
  X[, feats, drop = FALSE]
}
designs <- lapply(designs, function(z) { z$X <- align_X(z$X, all_feats); z })

if (length(all_feats) == 0) stop("No usable predictors after preprocessing.")

# ---------- Per-imputation CV over alpha grid ----------
set.seed(SEED)
alpha_results <- vector("list", length(ALPHA_GRID)); names(alpha_results) <- as.character(ALPHA_GRID)

for (a in ALPHA_GRID) {
  alpha_results[[as.character(a)]] <- lapply(seq_along(designs), function(k) {
    yk <- designs[[k]]$y
    Xk <- designs[[k]]$X
    if (length(yk) != nrow(Xk)) stop("y/X length mismatch in imputation ", k)
    
    cvfit <- glmnet::cv.glmnet(
      x = Xk, y = yk,
      family = "gaussian", alpha = a, nfolds = N_FOLDS,
      weights = rep(1/m, length(yk)), standardize = TRUE, parallel = FALSE
    )
    list(cvfit = cvfit,
         cvm_min = min(cvfit$cvm),
         lambda_min = cvfit$lambda.min,
         lambda_1se = cvfit$lambda.1se)
  })
}

alpha_perf <- tibble::tibble(
  alpha = ALPHA_GRID,
  avg_cvm = sapply(ALPHA_GRID, function(a) {
    mean(sapply(alpha_results[[as.character(a)]], `[[`, "cvm_min"))
  })
)
alpha_star <- alpha_perf$alpha[which.min(alpha_perf$avg_cvm)]
message(sprintf("Chosen alpha = %.2f (min average CV MSE across imputations).", alpha_star))

cvlist <- alpha_results[[as.character(alpha_star)]]

# ---------- Coefficients & selection frequencies (2-D matrix; robust) ----------
b1 <- coef(cvlist[[1]]$cvfit, s = cvlist[[1]]$lambda_1se)
feat_names <- rownames(b1)

coef_mat <- do.call(cbind, lapply(cvlist, function(obj) {
  b <- coef(obj$cvfit, s = obj$lambda_1se)
  as.numeric(b[feat_names, , drop = FALSE])
}))
rownames(coef_mat) <- feat_names
colnames(coef_mat) <- paste0("imp", seq_len(ncol(coef_mat)))

ix <- rownames(coef_mat) != "(Intercept)"
nonzero_freq <- rowMeans(coef_mat[ix, , drop = FALSE] != 0)
coef_median  <- apply(coef_mat[ix, , drop = FALSE], 1, median)
coef_mean    <- rowMeans(coef_mat[ix, , drop = FALSE])
coef_sd      <- apply(coef_mat[ix, , drop = FALSE], 1, sd)
coef_iqr     <- apply(coef_mat[ix, , drop = FALSE], 1, IQR)
coef_sign    <- sign(coef_mean)

summary_tbl <- tibble::tibble(
  feature     = rownames(coef_mat)[ix],
  sel_freq    = nonzero_freq,
  coef_median = coef_median,
  coef_mean   = coef_mean,
  coef_sd     = coef_sd,
  coef_iqr    = coef_iqr,
  sign_mean   = coef_sign
) %>% dplyr::arrange(dplyr::desc(sel_freq), dplyr::desc(abs(coef_median)))

readr::write_csv(alpha_perf, file.path(base_dir, "elasticnet_alpha_performance.csv"))
readr::write_csv(summary_tbl, file.path(base_dir, "elasticnet_selection_summary.csv"))
message("Saved alpha performance & selection summary.")

# ---------- Per-imputation metrics & predictions ----------
metrics_tbl <- purrr::imap_dfr(cvlist, function(obj, k) {
  Xk <- designs[[k]]$X; yk <- designs[[k]]$y
  yhat <- as.numeric(predict(obj$cvfit, newx = Xk, s = obj$lambda_1se))
  tibble::tibble(
    imp  = as.integer(k),
    n    = length(yk),
    mse  = mean((yk - yhat)^2),
    rmse = sqrt(mean((yk - yhat)^2)),
    r2   = 1 - sum((yk - yhat)^2) / sum((yk - mean(yk))^2)
  )
})
readr::write_csv(metrics_tbl, file.path(base_dir, "elasticnet_per_imputation_metrics.csv"))
message("Saved per-imputation metrics.")

# dump nonzero coefs per imputation
perimp_nonzero <- purrr::imap_dfr(cvlist, function(obj, k) {
  b <- as.matrix(coef(obj$cvfit, s = obj$lambda_1se))
  tibble::tibble(
    imp = as.integer(k),
    feature = rownames(b),
    coef = as.numeric(b)
  ) %>% dplyr::filter(feature != "(Intercept)", coef != 0)
})
readr::write_csv(perimp_nonzero, file.path(base_dir, "elasticnet_per_imp_nonzero_coefs.csv"))
message("Saved per-imputation nonzero coefficients.")

# predictions per imputation with ID if present
id_col <- if ("ID" %in% names(Dlong)) "ID" else NA_character_
pred_rows <- purrr::imap_dfr(cvlist, function(obj, k) {
  d <- designs[[k]]$data
  Xk <- designs[[k]]$X
  yhat <- as.numeric(predict(obj$cvfit, newx = Xk, s = obj$lambda_1se))
  tibble::tibble(
    imp = as.integer(k),
    row = seq_len(nrow(d)),
    ID  = if (!is.na(id_col) && id_col %in% names(d)) d[[id_col]] else NA,
    y   = d[[Y_POST]],
    yhat = yhat
  )
})
readr::write_csv(pred_rows, file.path(base_dir, "elasticnet_per_imp_predictions.csv"))
message("Saved per-imputation predictions.")

# ---------- Pooled OLS on stable support (Rubin's rules via miceadds::pool_mi) ----------
if (DO_STABLE_POOLS) {
  stable_feats <- summary_tbl %>% dplyr::filter(sel_freq >= STABILITY_THRESH) %>% dplyr::pull(feature)
  
  if (!length(stable_feats)) {
    warning("No features passed stability threshold (", STABILITY_THRESH, "). Skipping pooled OLS.")
  } else {
    message("Stable-support size: ", length(stable_feats), " (threshold = ", STABILITY_THRESH, ").")
    
    # Fit OLS on each imputation for SAME columns
    ols_list <- lapply(seq_along(designs), function(k) {
      d  <- designs[[k]]$data
      Xk <- designs[[k]]$X[, stable_feats, drop = FALSE]
      dk <- data.frame(y = d[[Y_POST]], as.data.frame(Xk), check.names = FALSE)
      fit <- lm(y ~ ., data = dk)
      list(coef = coef(fit), vcov = vcov(fit), names = names(coef(fit)))
    })
    
    # Align coefficient vectors/matrices across imputations
    all_coef_names <- Reduce(union, lapply(ols_list, function(z) z$names))
    
    Qlist <- lapply(ols_list, function(z) {
      b <- setNames(z$coef, z$names)
      b <- b[all_coef_names]
      b[is.na(b)] <- 0
      as.numeric(b)
    })
    
    Ulist <- lapply(ols_list, function(z) {
      V <- z$vcov
      M <- matrix(0, nrow = length(all_coef_names), ncol = length(all_coef_names),
                  dimnames = list(all_coef_names, all_coef_names))
      rn <- intersect(rownames(V), all_coef_names)
      cn <- intersect(colnames(V), all_coef_names)
      M[rn, cn] <- V[rn, cn]
      M
    })
    
    # Correct argument names: qhat=, u=
    pool_res <- miceadds::pool_mi(qhat = Qlist, u = Ulist)
    
    pooled_tbl <- tibble::tibble(
      term = all_coef_names,
      est  = pool_res$qbar,
      se   = sqrt(diag(pool_res$t)),
      tval = est / se,
      df   = pool_res$df,
      pval = 2 * pt(-abs(tval), df = df)
    ) %>% dplyr::arrange(pval)
    
    readr::write_csv(pooled_tbl, file.path(base_dir, "elasticnet_stable_support_pooled_OLS.csv"))
    message("Saved pooled OLS on stable support.")
  }
}

message("All done. Outputs in: ", normalizePath(base_dir, winslash = "/"))
