# === Numeric-only MICE; exclude genotype/Sex/Drug/Type/ID as predictors ===
need <- c("readr","dplyr","mice","naniar")
to_install <- setdiff(need, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(need, library, character.only = TRUE))

# 1) Load
path <- "E:/gwas project/casecase.csv"   # use forward slashes on Windows
stopifnot(file.exists(path))
df <- readr::read_csv(path, show_col_types = FALSE)

# 2) Identify columns
id_col   <- intersect("ID", names(df))
drug_col <- intersect("pharma.Drug", names(df))
type_col <- intersect("Type", names(df))

# all genotype.* (includes genotype.Sex and markers)
genotype_cols <- grep("^genotype\\.", names(df), value = TRUE)

# numeric blocks to impute
num_blocks <- c("^expressions\\.", "^concentrations\\.", "^cprs-r\\.",
                "^pharma\\.pre\\.", "^pharma\\.post\\.")
numeric_cols <- unique(unlist(lapply(num_blocks, function(p)
  grep(p, names(df), value = TRUE))))

# 3) Coerce ONLY numeric-block columns to numeric (preserves observed values)
numify <- function(x) suppressWarnings(as.numeric(x))
for (nm in numeric_cols) if (nm %in% names(df)) df[[nm]] <- numify(df[[nm]])

# 4) Quick audit BEFORE
cat("\n--- Missingness (numeric candidates) BEFORE ---\n")
print(naniar::miss_var_summary(df[, numeric_cols, drop = FALSE]) |>
        dplyr::arrange(dplyr::desc(pct_miss)) |> head(25))

# 5) Build mice specs
dat_imp <- df
meth <- make.method(dat_imp)
pred <- make.predictorMatrix(dat_imp)

# never impute non-numeric columns
non_numeric <- setdiff(names(dat_imp), numeric_cols)
meth[non_numeric] <- ""

# numeric columns: PMM only if there are at least 2 observed values
for (nm in numeric_cols) {
  obs <- sum(!is.na(dat_imp[[nm]]))
  meth[nm] <- if (obs >= 2) "pmm" else ""
}

# EXCLUDE these from predictor *and* target roles entirely
exclude_as_predictors <- unique(c(genotype_cols, drug_col, type_col, id_col))
if (length(exclude_as_predictors)) {
  pred[, exclude_as_predictors] <- 0  # do not use as predictors
  pred[exclude_as_predictors, ] <- 0  # and do not try to impute them
  meth[exclude_as_predictors] <- ""   # (redundant safety)
}

# also drop all-NA columns from prediction
all_na <- names(dat_imp)[colSums(!is.na(dat_imp)) == 0]
if (length(all_na)) { pred[, all_na] <- 0; pred[all_na, ] <- 0; meth[all_na] <- "" }

# no self-prediction
diag(pred) <- 0

# 6) Run MICE (observed cells never changed)
set.seed(20251105)
m <- 20; iters <- 20
imp <- mice(dat_imp, m = m, maxit = iters, method = meth,
            predictorMatrix = pred, printFlag = TRUE)

# 7) Save outputs
out_long <- "E:/gwas project/casecase_numeric_imputed_long_exclGD.csv"
readr::write_csv(complete(imp, action = "long", include = TRUE), out_long)
cat("\nSaved:", out_long, "\n")

out_1 <- "E:/gwas project/casecase_numeric_imputed_m1_exclGD.csv"
readr::write_csv(complete(imp, action = 1), out_1)
cat("Saved:", out_1, "\n")

# 8) Audit AFTER (on one completed set)
cat("\n--- Missingness (numeric candidates) AFTER (m=1) ---\n")
print(naniar::miss_var_summary(readr::read_csv(out_1, show_col_types = FALSE)[, numeric_cols, drop = FALSE]) |>
        dplyr::arrange(dplyr::desc(pct_miss)) |> head(25))
