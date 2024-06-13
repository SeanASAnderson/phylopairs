// This is a simple beta regression model
// The first column of the design matrix X is the intercept column and contains ones
// Each additional column contains observations for one predictor variable
// The model is g(mu.y) = XB, where g() is a link function, mu.y is the expected value of Ys given each observation of Xs...
// B is the vector of coefficients, and X is a design matrix

data {
  int<lower=0> N; // Number of observations (species pairs)
  int<lower=0> K; // Number of predictors, including the intercept
  vector<lower=0,upper=1>[N] Y; // Response variable, bounded between 0 and 1
  matrix[N, K] X; // Design matrix with intercept column and observations for each predictor variable
  real coef_mean; // mean for the normal prior on the coefficients
  real<lower=0> coef_sd; // sd for the normal prior on coefficients
  real<lower=0> phi_shape; // shape parameter for the gamma prior on phi
  real<lower=0> phi_rate; // rate parameter for the gamma prior on phi
  int<lower=1, upper=4> link_choice;  // 1 = logit, 2 = probit, 3 = cloglog, 4 = loglog
}

parameters {
  // note: while the beta distribution is defined by phi and mu, mu is calculated from the linear predictor in a beta regression
  // therefore, mu is defined in the model block and not in the parameters blcok 
  vector[K] Beta; // coefficient vector to be estimated
  real<lower=0> phi; // precision parameter for beta distribution
}

model {
  // Priors
  Beta ~ normal(0, 10); // Weakly informative priors for each coefficient
  phi ~ gamma(0.01, 0.01);

  // Calculate the linear predictor
  vector[N] eta = X*Beta;

  // Calculate the predicted values, mu, based on the link function
  vector[N] mu; // Mean of the Beta distribution
  if (link_choice == 1) {
    mu = inv_logit(eta);  // logit link
  } else if (link_choice == 2) {
    for (n in 1:N) {
      mu[n] = normal_cdf(eta[n], 0, 1);  // probit link
    }
  } else if (link_choice == 3) {
    mu = 1 - exp(-exp(eta));  // cloglog link
  } else if (link_choice == 4) {
    mu = exp(-exp(-eta));  // loglog link
  }
  
  // Likelihood
  Y ~ beta(mu * phi, (1 - mu) * phi);
}

generated quantities {
  vector[N] loglik = rep_vector(0.0, N); 
  vector[N] mu = inv_logit(X * Beta);

  // Calculate log likelihood in a vectorized manner
  for (n in 1:N) {
    loglik[n] = beta_lpdf(Y[n] | mu[n] * phi, (1 - mu[n]) * phi);
  }
}
