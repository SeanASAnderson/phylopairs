// This is a mixed model of linear regression
// In this model, there are fixed effects encoded in a design matrix X
// There are also pair-specific random effects (random intercepts)
// Pair-specific effects are not independent but covary according to the lineage-pair covariance matrix Cp
// The model, then, is Y = XB + u + e, u~N(0, sig2c*C) and e~N(0, sig2i*I)
// This is equivalent to Y = XB + eps where eps~N(0,(sig2c*C + sig2i*I)

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
  real sig2_loc; // location parameter of cauchy prior on sigma_resid
  real sig2_scale; // scale parameter of cauchy prior  on sigma_resid
}

parameters {
  vector[K] Beta; // coefficient vector to be estimated
  real<lower=0> sigma_resid; // standard deviation of residual errors
  real<lower=0> sig2_scale; // scaling parameter for variance in species-specific effects
  vector[N] pair_effects; // random effects (intercepts) for species pairs
}

model {
  // priors
  Beta ~ normal(coef_mean, coef_sd);
  sigma_resid ~ cauchy(sig2_loc, sig2_scale);
  sig2_scale ~ lognormal(sig2_mean, sig2_sd);
  pair_effects ~ multi_normal(rep_vector(0, N), sig2_scale * Cp);
  
  // Likelihood for response
  Y ~ normal(X*Beta + pair_effects, sigma_resid);
}

generated quantities {
  vector[N] y_pred;  // Predicted values of Y
  vector[N] loglik; // Log likelihood of each observation
  
  for (n in 1:N) {
    y_pred[n] = normal_rng(X[n] * Beta + pair_effects[n], sigma_resid);  // Generate predicted values
    loglik[n] = normal_lpdf(Y[n] | X[n]*Beta + pair_effects[n], sigma_resid);
  }
}
