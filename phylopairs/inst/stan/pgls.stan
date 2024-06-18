// This is a phylogenetic generalized least squares model
// The model is Y = XB + e, e~N(0, sig2_scale*C) where C is a lineage-pair covariance matrix (if analyzing species-pair data)
// or a phylogenetic covariance matrix (if analyzing species data)

data {
  int<lower=0> N; // Number of observations (lineage pairs)
  int<lower=0> K; // Number of predictors, including the intercept
  vector[N] Y; // Response variable
  matrix[N, K] X; // Design matrix with intercept column and observations for each predictor variable in subsequent columns
  matrix[N, N] Cp; // The lineage-pair covariance matrix
  real coef_mean; // mean for the normal prior on the coefficients
  real<lower=0> coef_sd; // sd for the normal prior on coefficients
  real sig2_mean; // mean of the lognormal prior on sig2_scale
  real sig2_sd; // sd of the lognormal prior prior on sig2_scale
}

parameters {
  vector[K] Coef; // coefficient vector to be estimated
  real<lower=0> sig2_scale; // scaling parameter for non-independent residuals
}

model {
  // priors
  Coef ~ normal(coef_mean, coef_sd);
  sig2_scale ~ lognormal(sig2_mean, sig2_sd);
  
  // Likelihood for response
  Y ~ multi_normal_cholesky(X * Coef, sqrt(sig2_scale) * cholesky_decompose(Cp));
}

generated quantities {
  vector[N] y_pred;  // Predicted values of Y
  vector[N] loglik;  // Log likelihood of each observation

  for (n in 1:N) {
    // Generate predicted value using the normal distribution
    y_pred[n] = normal_rng(X[n] * Coef, sqrt(sig2_scale * Cp[n, n]));

    // Calculate log likelihood for the observed value Y[n]
    loglik[n] = normal_lpdf(Y[n] | X[n] * Coef, sqrt(sig2_scale * Cp[n, n]));
  }
}
