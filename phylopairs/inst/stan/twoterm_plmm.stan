// This is a mixed model of linear regression fashioned after Castillo's model from (2007) Ecology and Evolution DOI: 10.1002/ece3.3093
// In this model, there are fixed effects encoded in a design matrix X (where the first column is the intercept column: a vector of 1s)
// There are also two classes of random effects, which are the species-specific effects and their covariance
// Species-specific effects are not independent but covary according to the phylogenetic covariance matrix C
// Each pair has two species-specific effect: one for "species 1" and one for "species 2"
// The model is Y = XB + Z1u1 + Z1u2 + e, where B is the vector of coefficients, X is a design matrix,
// u1 and u2 are vectors of species-specific random effects that covary according to C (u1 ~ MVN(0,C) and u2 ~ MVN(0,C))
// Z1 and Z2 are design matrices that map the two species-specific effects to each pair, and e is residual error
// Note that the random effects affect only the intercept: there are no random slopes.

data {
  int<lower=0> N; // Number of observations (lineage pairs)
  int<lower=0> M; // Number of species 1
  int<lower=0> P; // Number of species 2
  int<lower=0> K; // Number of predictors, including the intercept
  vector[N] Y; // Response variable
  matrix[N, K] X; // Design matrix with intercept column and observations for each predictor variable in subsequent columns
  matrix[M, M] C1; // The phylogenetic covariance matrix for species 1s
  matrix[P, P] C2; // The phylogenetic covariance matrix for species 2s
  matrix[N, M] Z1; // The design matrix for species-specific effects of 'species 1' in each pair
  matrix[N, P] Z2; // The design matrix for species-specific effects of 'species 2' in each pair
}

parameters {
  vector[K] Coef; // coefficient vector to be estimated
  real<lower=0> sigma_resid; // standard deviation of residual errors
  real<lower=0> sig2_scale1; // scaling parameter for variance in species-specific effects for species 1
  real<lower=0> sig2_scale2; // scaling parameter for variance in species-specific effects for species 2
  vector[M] sp1_effects; // random effects (intercepts) for 'species 1' in each pair
  vector[P] sp2_effects; // random effects (intercepts) for 'species 2' in each pair
  real coef_mean; // mean for the normal prior on the coefficients
  real<lower=0> coef_sd; // sd for the normal prior on coefficients
  real sig2_mean; // mean of the lognormal prior on sig2_scale
  real sig2_sd; // sd of the lognormal prior prior on sig2_scale
  real sig2_loc; // location parameter of cauchy prior on sigma_resid
  real sig2_sc; // scale parameter of cauchy prior  on sigma_resid
}

model {
  // priors
  //a ~ normal(0, 10);
  Coef ~ normal(0, 10);
  sigma_resid ~ cauchy(sig2_loc, sig2_sd);
  sig2_scale1 ~ lognormal(-1, 1);
  sig2_scale2 ~ lognormal(-1, 1);
  sp1_effects ~ multi_normal(rep_vector(0, M), sig2_scale1*C1);
  sp2_effects ~ multi_normal(rep_vector(0, P), sig2_scale2*C2);
  
  // Likelihood for response
  Y ~ normal(X*Coef + Z1*sp1_effects + Z2*sp2_effects, sigma_resid);
}

generated quantities {
  vector[N] loglik;  // Vector to store log-likelihood of each observation
  for (n in 1:N) {
    loglik[n] = normal_lpdf(Y[n] | (X[n] * Coef + Z1[n] * sp1_effects + Z2[n] * sp2_effects), sigma_resid);
  }
}


