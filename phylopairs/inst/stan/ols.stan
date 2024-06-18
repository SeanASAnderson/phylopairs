// This is a basic ols linear regression
// In this model, there are only fixed effects
// The model is Y = XB + e, where B is the vector of coefficients, X is a design matrix,
// and and e is residual error e~N(0, sigma_resid*I)

data {
  int<lower=0> N; // Number of observations (lineage pairs)
  int<lower=0> K; // Number of predictors, including the intercept
  vector[N] Y; // Response variable
  matrix[N, K] X; // Design matrix with intercept column and observations for each predictor variable in subsequent columns
  real coef_mean; // mean for the normal prior on the coefficients
  real<lower=0> coef_sd; // sd for the normal prior on coefficients
  real sig2_loc; // location parameter of cauchy prior on sigma_resid
  real sig2_sc; // scale parameter of cauchy prior  on sigma_resid
}

parameters {
  vector[K] Coef; // coefficient vector to be estimated
  real<lower=0> sigma_resid; // standard deviation of residual errors
}

model {
  // priors
  Coef ~ normal(coef_mean, coef_sd);
  sigma_resid ~ cauchy(sig2_loc, sig2_sc);

  // Likelihood for response
  Y ~ normal(X*Coef, sigma_resid);
}

generated quantities {
  vector[N] y_pred;  // Predicted values of Y
  vector[N] loglik; // Log likelihood of each observation
  
  for (n in 1:N) {
    y_pred[n] = normal_rng(X[n] * Coef, sigma_resid);  // Generate predicted values
    loglik[n] = normal_lpdf(Y[n] | X[n] * Coef, sigma_resid);  // Calculate log likelihood for each observation
  }
}

