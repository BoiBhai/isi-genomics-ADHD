## =========================
##      PLOTS ONLY
## =========================
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

# Defaults/helpers (safe if not already defined)
if (!exists("base_dir")) base_dir <- getwd()
if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)
if (!exists("PLOT_DPI")) PLOT_DPI <- 200
if (!exists("TOPN_PLOT")) TOPN_PLOT <- 30

safe_ggsave <- function(p, fn, width=8, height=5, dpi=PLOT_DPI) {
  fn <- file.path(base_dir, fn)
  tryCatch(
    ggsave(filename = fn, plot = p, width = width, height = height, dpi = dpi),
    error = function(e) message("Plot save failed: ", basename(fn), " -> ", e$message)
  )
}

## 1) Alpha performance (avg CV MSE) with chosen alpha highlighted
if (exists("alpha_perf") && exists("alpha_star")) {
  p_alpha <- ggplot(alpha_perf, aes(x = alpha, y = avg_cvm)) +
    geom_line() + geom_point(size = 2) +
    geom_point(data = dplyr::filter(alpha_perf, alpha == alpha_star), size = 3) +
    labs(title = "Average CV MSE vs alpha (across imputations)",
         x = expression(alpha), y = "Average CV MSE") +
    theme_minimal(base_size = 12)
  safe_ggsave(p_alpha, "plot_alpha_performance.png", width = 7, height = 4.5)
} else {
  message("Skip: alpha_perf/alpha_star not found.")
}

## 2) Selection frequency — top N features
if (exists("summary_tbl") && nrow(summary_tbl) > 0) {
  top_sel <- summary_tbl %>%
    arrange(desc(sel_freq), desc(abs(coef_median))) %>%
    slice_head(n = TOPN_PLOT)
  
  p_sel <- ggplot(top_sel, aes(x = reorder(feature, sel_freq), y = sel_freq)) +
    geom_col() +
    coord_flip() +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(title = paste0("Selection Frequency (Top ", nrow(top_sel), ")"),
         x = NULL, y = "Non-zero fraction across imputations") +
    theme_minimal(base_size = 11)
  safe_ggsave(p_sel, "plot_selection_frequency_top30.png", width = 8, height = 8)
} else {
  message("Skip: summary_tbl empty or missing (selection frequency).")
}

## 3) Coefficient median with IQR (25%–75%) — top N by |median|
if (exists("summary_tbl") && nrow(summary_tbl) > 0) {
  # Ensure q25/q75 exist; if not, approximate using IQR around median
  if (!all(c("q25","q75") %in% names(summary_tbl))) {
    summary_tbl <- summary_tbl %>%
      mutate(q25 = coef_median - coef_iqr/2,
             q75 = coef_median + coef_iqr/2)
  }
  top_coef <- summary_tbl %>%
    arrange(desc(abs(coef_median))) %>%
    slice_head(n = TOPN_PLOT)
  
  p_coef <- ggplot(top_coef,
                   aes(x = reorder(feature, abs(coef_median)),
                       y = coef_median)) +
    geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.2, alpha = 0.6) +
    geom_point(size = 2) +
    coord_flip() +
    labs(title = paste0("Coefficient Median with IQR (Top ", nrow(top_coef), ")"),
         x = NULL, y = "Median coefficient (IQR bars across imputations)") +
    theme_minimal(base_size = 11)
  safe_ggsave(p_coef, "plot_coef_median_IQR_top30.png", width = 8, height = 8)
} else {
  message("Skip: summary_tbl empty or missing (coef median/IQR).")
}

## 4) Overlay of CV curves (chosen alpha’s cvfits, all imputations)
if (exists("cvlist") && length(cvlist) > 0) {
  df_cv <- purrr::imap_dfr(cvlist, function(obj, k) {
    tibble(
      imp = as.integer(k),
      log_lambda = log(obj$cvfit$lambda),
      cvm = obj$cvfit$cvm
    )
  })
  df_min <- df_cv %>%
    group_by(imp) %>% slice_min(cvm, n = 1) %>% ungroup()
  
  p_cv <- ggplot(df_cv, aes(x = log_lambda, y = cvm, group = imp)) +
    geom_line(alpha = 0.25) +
    geom_point(data = df_min, size = 1.3, alpha = 0.6) +
    labs(title = "CV Curves per Imputation (chosen alpha)",
         x = "log(lambda)", y = "CV MSE") +
    theme_minimal(base_size = 12)
  safe_ggsave(p_cv, "plot_cv_curves_overlay.png", width = 7.5, height = 4.8)
} else {
  message("Skip: cvlist missing.")
}

## 5) Observed vs Predicted (all imputations stacked)
if (exists("pred_rows") && all(c("y","yhat") %in% names(pred_rows))) {
  p_scatter <- ggplot(pred_rows, aes(x = y, y = yhat)) +
    geom_point(alpha = 0.25) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = "Observed vs Predicted (all imputations)",
         x = "Observed", y = "Predicted") +
    theme_minimal(base_size = 12)
  safe_ggsave(p_scatter, "plot_obs_vs_pred_allimps.png", width = 6.5, height = 6)
} else {
  message("Skip: pred_rows missing or lacks y/yhat.")
}

## 6) Residual density (all imputations)
if (exists("pred_rows") && all(c("y","yhat") %in% names(pred_rows))) {
  pred_rows$resid <- pred_rows$y - pred_rows$yhat
  p_res <- ggplot(pred_rows, aes(x = resid)) +
    geom_density() +
    labs(title = "Residual Density (all imputations)", x = "Residual", y = "Density") +
    theme_minimal(base_size = 12)
  safe_ggsave(p_res, "plot_residual_density_allimps.png", width = 7, height = 4.5)
} else {
  message("Skip: pred_rows missing or lacks y/yhat (residuals).")
}

## 7) Coefficient paths (first imputation’s glmnet.fit)
if (exists("cvlist") && length(cvlist) >= 1 &&
    !is.null(cvlist[[1]]$cvfit$glmnet.fit)) {
  png(file.path(base_dir, "plot_glmnet_paths_imp1.png"),
      width = 1200, height = 800, res = PLOT_DPI)
  plot(cvlist[[1]]$cvfit$glmnet.fit, xvar = "lambda", label = FALSE)
  abline(v = log(cvlist[[1]]$lambda_1se), lty = 2, col = "grey40")
  dev.off()
} else {
  message("Skip: glmnet path plot not available (cvlist[[1]]$cvfit$glmnet.fit missing).")
}
