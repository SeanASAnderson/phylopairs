// This is a mixed model version of beta regression
// In this model, there are fixed effects encoded in a design matrix X
// The first column of X is the intercept column and contains ones
// Each additional column contains observations for one predictor variable
// There is also a class of random effect, which are random intercepts for each pair
// Pair-specific effects are not independent but covary according to the lineage-pair covariance matrix Cp 
// (or to the phylogenetic covariance matrix C if using species data)
// The model, then, is g(mu.y) = XB + u, where g() is a link function, mu.y is the expected value of Ys given each observation of Xs...
// B is the vector of coefficients, X is a design matrix, and u is the vector of pair-specific random effects ("pair_effects")

data {
  int<lower=0> N; // Number of observations (lineage pairs)
  int<lower=0> K; // Number of predictors, including the intercept
  vector<lower=0,upper=1>[N] Y; // Response variable, bounded between 0 and 1
  matrix[N, K] X; // Design matrix with intercept column and observations for each predictor variable in subsequent columns
  matrix[N, N] Cp; // The lineage-pair covariance matrix
  real coef_mean; // mean for the normal prior on the coefficients
  real<lower=0> coef_sd; // sd for the normal prior on coefficients
  real<lower=0> phi_shape; // shape parameter for the gamma prior on phi
  real<lower=0> phi_rate; // rate parameter for the gamma prior on phi
  real sig2_mean; // mean of the lognormal prior on sig2_scale
  real sig2_sd; // sd of the lognormal prior prior on sig2_scale
  int<lower=1, upper=4> link_choice;  // 1 = logit, 2 = probit, 3 = cloglog, 4 = loglog
}

parameters {
  // note: this parameterization of beta distribution is defined by phi and mu
  // however, mu is defined in the model block, not the parameters block
  // this is because mu is calculated from the linear predictor
  // thus, we estimate the parameters for the linear predictor that give you mu, not mu itself
  vector[K] Coef; // coefficient vector to be estimated
  real<lower=0> phi; // precision parameter for beta distribution
  real<lower=0> sig2_scale; // scaling parameter for variance in species-specific effects
  vector[N] pair_effects; // random effects (intercepts) for species pairs
}

model {
  // priors
  Coef ~ normal(coef_mean, coef_sd);
  phi ~ gamma(phi_shape, phi_rate);
  sig2_scale ~ lognormal(sig2_mean, sig2_sd);
  pair_effects ~ multi_normal(rep_vector(0, N), sig2_scale*Cp);

  // Calculate linear predictor
  vector[N] mu;
  vector[N] eta = X*Coef + pair_effects; 

  // Apply the chosen link function to get the predicted value of each y
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

  // Likelihood for response
  Y ~ beta(mu * phi, (1 - mu) * phi);
}

generated quantities {
  vector[N] loglik; 
  vector[N] mu; // Mean of the Beta distribution

  // Calculate the predicted values based on the link function
  if (link_choice == 1) {
    mu = inv_logit(X * Coef + pair_effects);  // logit link
  } else if (link_choice == 2) {
    for (n in 1:N) {
      mu[n] = normal_cdf(X[n] * Coef + pair_effects[n], 0, 1);  // probit link
    }
  } else if (link_choice == 3) {
    mu = 1 - exp(-exp(X * Coef + pair_effects));  // cloglog link
  } else if (link_choice == 4) {
    mu = exp(-exp(-(X * Coef + pair_effects)));  // loglog link
  }

  // Calculate log likelihood in a vectorized manner
  for (n in 1:N) {
    loglik[n] = beta_lpdf(Y[n] | mu[n] * phi, (1 - mu[n]) * phi);
  }
}

