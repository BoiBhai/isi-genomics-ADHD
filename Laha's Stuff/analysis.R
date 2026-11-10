#!/usr/bin/env Rscript
library(tidyverse)
library(broom)
library(mice)
library(sandwich)
library(lmtest)
library(glmnet)

input_file <- "./cleaned_joined_data [2025-11-05T202120.986+0530].csv"
data <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(across(where(is.character), ~ as_factor(.x)))

OUT_DIR <- "./out"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, showWarnings = FALSE)

## Minor allele frequency (MAF) calculation
# Identify genotype columns (those that start with "genotype.")
geno_metrics <- data %>%
  select(starts_with("genotype.M")) %>%
  reframe(across(
    everything(),
    ~ {
      # Extract alleles from genotype strings
      alleles <- unlist(str_split(na.omit(.x), ""))
      if (length(alleles) == 0) {
        return(NA_real_) # No data to calculate MAF
      }
      allele_counts <- table(alleles)
      total_alleles <- sum(allele_counts)
      freqs <- allele_counts / total_alleles
      if (length(freqs) < 2) {
        return(0) # Monomorphic site
      }
      major_allele <- names(freqs)[which.max(freqs)]
      minor_allele <- names(freqs)[which.min(freqs)]
      maf <- freqs[minor_allele]

      genotype_counts <- table(.x)
      dominant_count <- as.numeric(genotype_counts[paste0(major_allele, major_allele)])
      dominant_count <- ifelse(is.na(dominant_count), 0, dominant_count)
      hetero_count1 <- as.numeric(genotype_counts[paste0(major_allele, minor_allele)])
      hetero_count1 <- ifelse(is.na(hetero_count1), 0, hetero_count1)
      hetero_count2 <- as.numeric(genotype_counts[paste0(minor_allele, major_allele)])
      hetero_count2 <- ifelse(is.na(hetero_count2), 0, hetero_count2)
      hetero_count <- hetero_count1 + hetero_count2
      recessive_count <- as.numeric(genotype_counts[paste0(minor_allele, minor_allele)])
      recessive_count <- ifelse(is.na(recessive_count), 0, recessive_count)

      # Test for Hardy-Weinberg Equilibrium (HWE)
      n <- sum(genotype_counts)
      hwe_p_value <- chisq.test(c(
        dominant_count, hetero_count, recessive_count
      ), p = c(
        (1 - maf)^2,
        2 * maf * (1 - maf),
        maf^2
      ))$p.value

      tibble(value = c(
        maf, major_allele, minor_allele,
        dominant_count / n, hetero_count / n, recessive_count / n,
        hwe_p_value
      ))
    }
  )) %>%
  t() %>%
  as_tibble(rownames = "Marker") %>%
  set_names(c(
    "Marker", "MAF", "Major_Allele", "Minor_Allele",
    "Dominant_Freq", "Heterozygous_Freq", "Recessive_Freq",
    "HWE_P_Value"
  )) %>%
  mutate(
    Marker = str_remove(Marker, "genotype."),
    MAF = as.numeric(MAF),
    Dominant_Freq = as.numeric(Dominant_Freq),
    Heterozygous_Freq = as.numeric(Heterozygous_Freq),
    Recessive_Freq = as.numeric(Recessive_Freq),
    HWE_P_Value = as.numeric(HWE_P_Value)
  ) %>%
  mutate(
    HWE_Holds = ifelse(HWE_P_Value < 0.05, "No", "Yes")
  ) %>%
  arrange(MAF) %>%
  write_csv(file.path(OUT_DIR, "genotype_metrics.csv"))

## Exclude columns where HWE does not hold
data <- data %>%
  select(
    -any_of(
      geno_metrics %>%
        filter(HWE_Holds == "No") %>%
        pull(Marker) %>%
        paste0("genotype.", .)
    )
  )

## Ordinalise genotype columns based on major and minor alleles
data <- data %>%
  mutate(across(
    .cols = starts_with("genotype.M"),
    .fns = ~ {
      marker <- str_remove(cur_column(), "^genotypes?\\.")
      mm <- geno_metrics %>% filter(Marker == marker)
      if (nrow(mm) == 0) {
        return(rep(NA, length(.x)))
      }

      minor <- mm$Minor_Allele
      major <- mm$Major_Allele

      genotypes <- c(
        paste0(major, major), paste0(major, minor), paste0(minor, minor)
      )
      factor(.x, levels = genotypes, ordered = TRUE)
    },
    .names = "{.col}"
  ))

## Summary statistics for numeric columns
summary_stats <- data %>%
  select(where(is.numeric), Type) %>%
  pivot_longer(
    cols = where(is.numeric),
    names_to = "Variable", values_to = "Value"
  ) %>%
  group_by(Variable, Type) %>%
  summarise(
    n = sum(!is.na(Value)),
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    median = median(Value, na.rm = TRUE),
    siqr = IQR(Value, na.rm = TRUE),
    min = min(Value, na.rm = TRUE),
    max = max(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n > 0) %>%
  arrange(Variable, Type) %>%
  write_csv(file.path(OUT_DIR, "summary_statistics.csv"))

## Do case-control t-tests for expressions and concentrations
t_test_results <- data %>%
  select(starts_with("expressions"), starts_with("concentrations"), Type) %>%
  pivot_longer(
    cols = where(is.numeric),
    names_to = "Variable", values_to = "Value"
  ) %>%
  group_by(Variable) %>%
  summarise(
    test = tidy(t.test(
      Value ~ Type,
      data = pick(Value, Type),
    )),
    .groups = "drop"
  ) %>%
  unnest_wider(test) %>%
  mutate(
    d = estimate,
    mean1 = estimate1,
    mean2 = estimate2
  ) %>%
  mutate(
    p.adjusted.BH = p.adjust(p.value, method = "BH"),
    p.adjusted.Holm = p.adjust(p.value, method = "holm")
  ) %>%
  select(Variable, mean1, mean2, d, statistic, p.value, p.adjusted.BH, p.adjusted.Holm) %>%
  arrange(p.adjusted.BH) %>%
  write_csv(file.path(OUT_DIR, "t_test_results.csv"))

## Make a volcano plot for t-test results
volcano_data <- t_test_results %>%
  mutate(
    neg_log10_p = -log10(p.value),
    signif = case_when(
      p.adjusted.BH < 0.05 & abs(d) >= 1 ~ "Significant",
      TRUE ~ "Not Significant"
    )
  )
# create volcano plot
volcano_plot <- ggplot(volcano_data, aes(x = d, y = neg_log10_p, color = signif)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey50")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  labs(
    x = "Effect size (difference in means)",
    y = expression(-log[10](p.value)),
    title = "Volcano plot of t-test results",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  # annotate significant points to avoid overplotting; requires ggrepel package
  ggrepel::geom_text_repel(
    data = dplyr::filter(volcano_data, p.adjusted.BH < 0.05 & abs(d) >= 1),
    aes(label = Variable),
    size = 3,
    max.overlaps = 10,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.size = 0.2
  )

# save plot files
ggsave(file.path(OUT_DIR, "volcano_plot.pdf"), volcano_plot, bg = "transparent")

## Implement a MANOVA on expression and concentration data
manova_data <- data %>%
  select(starts_with("expressions"), starts_with("concentrations"), Type)

# Impute separately within each Type and return a single completed dataset by averaging imputations
manova_data_imputed <- manova_data %>%
  mutate(.row_id = row_number()) %>%
  group_by(Type) %>%
  group_split() %>%
  purrr::map_dfr(function(df_grp) {
    row_ids <- df_grp$.row_id
    df_work <- df_grp %>% select(-.row_id, -Type)

    if (!any(is.na(df_work))) {
      # no missing values -> return as-is
      return(bind_cols(df_work, Type = df_grp$Type, .row_id = row_ids))
    }

    imp <- mice::mice(df_work, m = 5, printFlag = FALSE, seed = 2402)
    long <- mice::complete(imp, action = "long", include = FALSE) %>% select(-.imp)

    averaged <- long %>%
      group_by(.id) %>%
      summarise(across(everything(), ~ if (is.numeric(.x)) mean(.x) else dplyr::first(.x)), .groups = "drop") %>%
      arrange(.id) %>%
      select(-.id)

    bind_cols(averaged, Type = df_grp$Type, .row_id = row_ids)
  }) %>%
  arrange(.row_id) %>%
  select(-.row_id)

manova_imputed <- manova(
  as.matrix(manova_data_imputed %>% select(-Type)) ~ Type,
  data = manova_data_imputed
)
capture.output(print(manova_imputed),
  file = file.path(OUT_DIR, "manova_imputed_details.txt")
)
capture.output(summary(manova_imputed),
  file = file.path(OUT_DIR, "manova_imputed_summary.txt")
)

regression_data <- data %>%
  mutate(
    Sex = factor(case_when(
      genotype.Sex == "M" ~ "Male",
      genotype.Sex == "F" ~ "Female",
      TRUE ~ "Unknown"
    ))
  ) %>%
  mutate(
    Sex = relevel(Sex, ref = "Unknown")
  )

## Trait modelling with HC3 sandwich estimator for robust standard errors

# First, look at the column names for expressions, concentrations, and cprs-r
expression_columns <- names(data) %>%
  str_subset("^expressions\\.")
names(expression_columns) <- expression_columns
concentration_columns <- names(data) %>%
  str_subset("^concentrations\\.")
names(concentration_columns) <- concentration_columns
cprs_columns <- names(data) %>%
  str_subset("^cprs-r\\.")
names(cprs_columns) <- cprs_columns


lm_hc3 <- function(formula, data) {
  model <- lm(formula, data = data)
  robust_fit <- lmtest::coeftest(
    model,
    vcov = sandwich::vcovHC(model, type = "HC3")
  )
  tidy(robust_fit)
}

# Traits vs Concentrations
purrr::map(cprs_columns, function(trait) {
  formula <- as.formula(
    paste(
      paste0("`", trait, "`"),
      "~ Sex + ",
      paste(
        concentration_columns |> map_chr(~ paste0("`", .x, "`")),
        collapse = " + "
      )
    )
  )
  output_file <- file.path(
    OUT_DIR, paste0("robust_regression_", trait, "_vs_Concentrations.csv")
  )
  lm_hc3(formula, regression_data) %>% write_csv(output_file)
})

# Concentrations vs Expressions
purrr::map(concentration_columns, function(chem) {
  formula <- as.formula(
    paste(
      paste0("`", chem, "`"),
      "~ Sex + ",
      paste(
        expression_columns |> map_chr(~ paste0("`", .x, "`")),
        collapse = " + "
      )
    )
  )
  output_file <- file.path(
    OUT_DIR, paste0("robust_regression_", chem, "_vs_Expressions.csv")
  )
  lm_hc3(formula, regression_data) %>% write_csv(output_file)
})

## cross-validated prediction of Type using concentrations using cv.glmnet
cv_data <- regression_data %>%
  select(all_of(concentration_columns), Type) %>%
  drop_na()

# Fit cross-validated elastic net model
cv_glmnet_model <- cv.glmnet(
  x = {
    cv_data %>%
      select(all_of(concentration_columns)) %>%
      as.matrix()
  },
  y = {
    (
      cv_data %>%
        pull(Type) == "Case"
    ) %>%
      as.numeric()
  },
  family = "binomial",
  alpha = 0.5, # Elastic net (0.5 = equal mix of L1 and L2)
  nfolds = 5,
  seed = 2402
)

# Extract results and save
cv_results <- tibble(
  lambda = cv_glmnet_model$lambda,
  cvm = cv_glmnet_model$cvm,
  cvsd = cv_glmnet_model$cvsd,
  nzero = cv_glmnet_model$nzero
)

write_csv(cv_results, file.path(OUT_DIR, "cv_glmnet_lambda_results_concentrations.csv"))

# Plot cross-validation error
cv_plot <- plot(cv_glmnet_model)
ggsave(file.path(OUT_DIR, "cv_glmnet_plot.pdf"), cv_plot, bg = "transparent")

# Extract coefficients at lambda.min
coefficients_min <- coef(cv_glmnet_model, s = "lambda.min") %>%
  as.matrix() %>%
  as_tibble(rownames = "Variable") %>%
  set_names(c("Variable", "Coefficient")) %>%
  filter(Variable != "(Intercept)") %>%
  arrange(desc(abs(Coefficient))) %>%
  write_csv(file.path(OUT_DIR, "cv_glmnet_coefficients_concentrations.csv"))

## cross-validated prediction of Type using expressions using cv.glmnet
cv_data <- regression_data %>%
  select(all_of(expression_columns), Type) %>%
  drop_na()

# Fit cross-validated elastic net model
cv_glmnet_model <- cv.glmnet(
  x = {
    cv_data %>%
      select(all_of(expression_columns)) %>%
      as.matrix()
  },
  y = {
    (
      cv_data %>%
        pull(Type) == "Case"
    ) %>% as.numeric()
  },
  family = "binomial",
  alpha = 0.5, # Elastic net (0.5 = equal mix of L1 and L2)
  nfolds = 5,
  seed = 2402
)

# Extract results and save
cv_results <- tibble(
  lambda = cv_glmnet_model$lambda,
  cvm = cv_glmnet_model$cvm,
  cvsd = cv_glmnet_model$cvsd,
  nzero = cv_glmnet_model$nzero
)

write_csv(cv_results, file.path(OUT_DIR, "cv_glmnet_lambda_results_expressions.csv"))
# Plot cross-validation error
cv_plot <- plot(cv_glmnet_model)
ggsave(file.path(OUT_DIR, "cv_glmnet_plot_expressions.pdf"), cv_plot, bg = "transparent")
# Extract coefficients at lambda.min
coefficients_min <- coef(cv_glmnet_model, s = "lambda.min") %>%
  as.matrix() %>%
  as_tibble(rownames = "Variable") %>%
  set_names(c("Variable", "Coefficient")) %>%
  filter(Variable != "(Intercept)") %>%
  arrange(desc(abs(Coefficient))) %>%
  write_csv(file.path(OUT_DIR, "cv_glmnet_coefficients_expressions.csv"))

## Multi-stage Bayesian modelling using Stan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Prepare data for Stan model
stan_data_raw <- regression_data %>%
  select(all_of(c(expression_columns, concentration_columns, cprs_columns)))

# Create missing data indicators (1 if observed, 0 if missing)
concentrations_matrix <- stan_data_raw %>%
  select(all_of(concentration_columns)) %>%
  as.matrix()
concentrations_ins <- (!is.na(concentrations_matrix)) * 1L

expressions_matrix <- stan_data_raw %>%
  select(all_of(expression_columns)) %>%
  as.matrix()
expressions_ins <- (!is.na(expressions_matrix)) * 1L

traits_matrix <- stan_data_raw %>%
  select(all_of(cprs_columns)) %>%
  as.matrix()

# Replace NAs with 0 for Stan (will be imputed in the model)
concentrations_matrix[is.na(concentrations_matrix)] <- 0
expressions_matrix[is.na(expressions_matrix)] <- 0
traits_matrix[is.na(traits_matrix)] <- 0

# Prepare Stan data list
stan_data <- list(
  N = nrow(stan_data_raw),
  T = length(cprs_columns),
  C = length(concentration_columns),
  E = length(expression_columns),
  traits = traits_matrix,
  concentrations = concentrations_matrix,
  concentrations_ins = concentrations_ins,
  expressions = expressions_matrix,
  expressions_ins = expressions_ins
)

cat("\nStan Data Summary:\n")
cat("  Number of observations (N):", stan_data$N, "\n")
cat("  Number of traits (T):", stan_data$T, "\n")
cat("  Number of concentrations (C):", stan_data$C, "\n")
cat("  Number of gene expressions (E):", stan_data$E, "\n")
cat("  Missing concentrations:", sum(stan_data$concentrations_ins == 0), "\n")
cat("  Missing expressions:", sum(stan_data$expressions_ins == 0), "\n")

# Compile Stan model
cat("\nCompiling Stan model...\n")
stan_model <- stan_model(file = "./multi-stage.stan")

# Fit the model with reduced iterations for stability/speed
cat("Fitting multi-stage Stan model...\n")
fit <- sampling(
  stan_model,
  data = stan_data,
  iter = 2000,
  warmup = 1000,
  chains = 5,
  cores = parallel::detectCores(),
  seed = 2402
)

# Extract and save results
cat("\nExtracting model results...\n")

# Summary statistics
stan_summary <- fit %>%
  summary() %>%
  .$summary %>%
  as_tibble(rownames = "Parameter") %>%
  write_csv(file.path(OUT_DIR, "stan_model_summary.csv"))

# Extract imputed concentrations
concentrations_imputed <- rstan::extract(fit, pars = "concentrations_imp")[[1]]
concentration_imp_mean <- apply(concentrations_imputed, 2:3, mean)
concentration_imp_sd <- apply(concentrations_imputed, 2:3, sd)

concentration_imp_results <- tibble(
  Observation = rep(seq_len(nrow(concentration_imp_mean)), times = ncol(concentration_imp_mean)),
  Concentration = rep(concentration_columns, each = nrow(concentration_imp_mean)),
  Mean = as.vector(concentration_imp_mean),
  SD = as.vector(concentration_imp_sd)
) %>%
  write_csv(file.path(OUT_DIR, "stan_imputed_concentrations.csv"))

# Extract imputed expressions
expressions_imputed <- rstan::extract(fit, pars = "expressions_imp")[[1]]
expressions_imp_mean <- apply(expressions_imputed, 2:3, mean)
expressions_imp_sd <- apply(expressions_imputed, 2:3, sd)

expressions_imp_results <- tibble(
  Observation = rep(seq_len(nrow(expressions_imp_mean)), times = ncol(expressions_imp_mean)),
  Expression = rep(expression_columns, each = nrow(expressions_imp_mean)),
  Mean = as.vector(expressions_imp_mean),
  SD = as.vector(expressions_imp_sd)
) %>%
  write_csv(file.path(OUT_DIR, "stan_imputed_expressions.csv"))

# Extract regression coefficients for concentrations -> traits
beta_conc <- rstan::extract(fit, pars = "beta_conc")[[1]]
beta_conc_mean <- apply(beta_conc, 2:3, mean)
beta_conc_results <- tibble(
  Concentration = rep(concentration_columns, times = length(cprs_columns)),
  Trait = rep(cprs_columns, each = length(concentration_columns)),
  Coefficient = as.vector(t(beta_conc_mean))
) %>%
  write_csv(file.path(OUT_DIR, "stan_beta_concentrations_traits.csv"))

# Extract regression coefficients for expressions -> concentrations
beta_expr <- rstan::extract(fit, pars = "beta_expr")[[1]]
beta_expr_mean <- apply(beta_expr, 2:3, mean)
beta_expr_results <- tibble(
  Expression = rep(expression_columns, times = length(concentration_columns)),
  Concentration = rep(concentration_columns, each = length(expression_columns)),
  Coefficient = as.vector(t(beta_expr_mean))
) %>%
  write_csv(file.path(OUT_DIR, "stan_beta_expressions_concentrations.csv"))

# Extract variance components
sigma_traits <- rstan::extract(fit, pars = "sigma_traits")[[1]]
sigma_concentrations <- rstan::extract(fit, pars = "sigma_concentrations")[[1]]

sigma_traits_mean <- apply(sigma_traits, 2:3, mean)
sigma_conc_mean <- apply(sigma_concentrations, 2:3, mean)

sigma_results <- tibble(
  Variable = c(
    rep(cprs_columns, each = nrow(sigma_traits_mean)),
    rep(concentration_columns, each = nrow(sigma_conc_mean))
  ),
  Observation = c(
    rep(seq_len(nrow(sigma_traits_mean)), times = length(cprs_columns)),
    rep(seq_len(nrow(sigma_conc_mean)), times = length(concentration_columns))
  ),
  Type = c(
    rep("Trait", length(as.vector(sigma_traits_mean))),
    rep("Concentration", length(as.vector(sigma_conc_mean)))
  ),
  Sigma = c(as.vector(sigma_traits_mean), as.vector(sigma_conc_mean))
) %>%
  write_csv(file.path(OUT_DIR, "stan_variance_components.csv"))

# Extract correlation matrices
Rho_traits <- rstan::extract(fit, pars = "Rho_traits")[[1]]
Rho_concentrations <- rstan::extract(fit, pars = "Rho_concentrations")[[1]]

Rho_traits_mean <- apply(Rho_traits, 2:3, mean)
Rho_conc_mean <- apply(Rho_concentrations, 2:3, mean)

# Save correlation matrices as CSV
write.csv(
  Rho_traits_mean,
  file.path(OUT_DIR, "stan_correlation_matrix_traits.csv"),
  row.names = cprs_columns
)
write.csv(
  Rho_conc_mean,
  file.path(OUT_DIR, "stan_correlation_matrix_concentrations.csv"),
  row.names = concentration_columns
)

# Save Stan fit object for diagnostics
saveRDS(fit, file.path(OUT_DIR, "stan_fit_object.rds"))

cat("Stan model fitting completed successfully!\n")
cat("Results saved to:", OUT_DIR, "\n")
