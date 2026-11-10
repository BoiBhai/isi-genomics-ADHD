data {
  int<lower=0> N; // number of observations
  int<lower=0> T; // number of traits
  int<lower=0> C; // number of concentrations
  int<lower=0> E; // number of gene expressions
  matrix[N, T] traits; // matrix of observed traits
  matrix[N, C] concentrations; // matrix of (observed) concentrations
  array[N, C] int<lower=0, upper=1> concentrations_ins; // array of indicators for missing concentrations (1 if observed, 0 if missing)
  matrix[N, E] expressions; // matrix of observed gene expressions
  array[N, E] int<lower=0, upper=1> expressions_ins; // array of indicators for missing expressions (1 if observed, 0 if missing)
}
parameters {
  matrix[N, C] concentrations_imp; // matrix of imputed concentrations
  matrix[N, E] expressions_imp; // matrix of imputed gene expressions
  matrix[C, T] beta_conc; // regression coefficients for concentrations (mean)
  matrix[E, C] beta_expr; // regression coefficients for expressions (mean)
  matrix[C, T] gamma_conc; // regression coefficients for concentrations (std dev)
  matrix[E, C] gamma_expr; // regression coefficients for expressions (std dev)
  
  // Hierarchical variance components
  corr_matrix[T] Rho_traits; // correlation matrix of the traits
  array[N] vector<lower=0>[T] sigma_traits; // standard deviations of the traits
  corr_matrix[C] Rho_concentrations; // correlation matrix of the concentrations
  array[N] vector<lower=0>[C] sigma_concentrations; // standard deviations of the concentrations
  
  cov_matrix[T] W_traits; // covariance matrix of the log standard deviations of the traits
  cov_matrix[C] W_concentrations; // covariance matrix of the log standard deviations of the concentrations
}
transformed parameters {
  matrix[N, C] concentrations_full;
  matrix[N, E] expressions_full;
  
  concentrations_full = concentrations .* to_matrix(concentrations_ins)
                        + concentrations_imp
                          .* (1 - to_matrix(concentrations_ins));
  expressions_full = expressions .* to_matrix(expressions_ins)
                     + expressions_imp .* (1 - to_matrix(expressions_ins));
  
  array[N] cov_matrix[T] Sigma_traits;
  for (i in 1 : N) {
    Sigma_traits[i] = Rho_traits .* (sigma_traits[i] * sigma_traits[i]');
  }
  array[N] cov_matrix[C] Sigma_concentrations;
  for (i in 1 : N) {
    Sigma_concentrations[i] = Rho_concentrations
                              .* (sigma_concentrations[i]
                                  * sigma_concentrations[i]');
  }
}
model {
  for (i in 1 : N) {
    traits[i] ~ multi_normal(concentrations_full[i] * beta_conc,
                             Sigma_traits[i]);
  }
  
  for (i in 1 : N) {
    concentrations_full[i] ~ multi_normal(expressions_full[i] * beta_expr,
                                          Sigma_concentrations[i]);
  }
  
  for (t in 1 : T) {
    beta_conc[ : , t] ~ normal(0, 5);
    gamma_conc[ : , t] ~ normal(0, 2);
  }
  
  Rho_traits ~ lkj_corr(2.0);
  for (i in 1 : N) {
    log(sigma_traits[i]) ~ multi_normal(concentrations_full[i] * gamma_conc,
                                        W_traits);
  }
  Rho_concentrations ~ lkj_corr(2.0);
  for (i in 1 : N) {
    log(sigma_concentrations[i]) ~ multi_normal(
                                                expressions_full[i]
                                                * gamma_expr,
                                                W_concentrations);
  }
  W_traits ~ inv_wishart(T + 1, diag_matrix(rep_vector(1, T)));
  W_concentrations ~ inv_wishart(C + 1, diag_matrix(rep_vector(1, C)));
}
