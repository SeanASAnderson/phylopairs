#' betareg_stan

#' @description Fits one of two beta regression models to a dataset in a Bayesian 
#'   framework using Stan software via the 'rstan' package. Response variable must be
#'   bounded between 0 and 1. Users can fit either (1) a basic beta regression model or a 
#'   (2) beta mixed-effects model in which there are covarying random intercepts in the 
#'   linear predictor. For the latter, users must supply a covariance matrix (the 
#'   simpler beta regression model is otherwise run). In both models, users can choose 
#'   one of four link functions. See details. 

#' @details
#'   The two models that can be chosen differ in terms of their linear predictor 
#'   on the scale of the link function. These predictors are
#'   1. Simple beta regression: \code{g(u) = XB}
#'   2. Beta mixed-effects regression: \code{g(u) = XB + u, u~N(0, physig*C)}, where C is a 
#'      lineage-pair covariance matrix (if analyzing lineage-pair data) or a 
#'      phylogenetic vcv matrix (if analyzing species data) and physig is the sigma^2 scaling 
#'      parameter for C. 
#'   
#'   NOTE: For model 2, the user must supply a covariance matrix.
#'   # Prior Distributions for Model Parameters #
#'   The underlying stan models assume the following prior distributions for model parameters
#'   1. Regression coefficients: Gaussian prior (users can set prior mean and sd)
#'   2. Precision parameter `phi`: gamma prior (users can set prior shape and rate of gamma)
#'   3. `physig`: lognormal prior (users can set prior mean and sd)

#' @param des A vector of predictor variable observations OR, in the case of multiple predictors, a matrix in which each column is a vector of observations of a given predictor. 'betareg.stan' adds a column of 1s to make this a design matrix whose first column corresponds to the model intercept (unless such a column already exists).
#' @param y A vector of response variable observations. 
#' @param cov Covariance matrix for model residuals (a Cp matrix if analyzing lineage-pair data or a phylogenetic vcv matrix if analyzing bounded species data).
#' @param link The link function to be used. Default is "logit". Other possible choices are "probit", "cloglog" (complementary log-log), and "loglog".
#' @param iter Number of iterations to run on each chain; defaults to 6000.
#' @param chains Number of chains to run; defaults to 4.
#' @param cores Number of cores to be used; defaults to 4 (one chain per core).
#' @param coef.u Mean of the Gaussian prior for each preditor variable coefficient; defaults to 0.
#' @param coef.sd SD of the Gaussian prior for each preditor variable coefficient; defaults to 10.
#' @param phi.shape Shape parameter for gamma prior of beta distribution's precision parameter (phi); defaults to 0.01.
#' @param phi.rate Rate parameter for gamma prior of beta distribution's precision parameter (phi); defaults to 0.01.
#' @param physig2.u Mean of the lognormal prior for the scale of the residual covariance; defaults to -1.
#' @param physig2.sd SD of the lognormal prior for the scale of the residual covariance; defaults to 1.

#' @return A list containing two elements: (1) the posterior distribution of beta parameters, and (2) the log-likelihood of the posteriors for use in downstream analyses (e.g. the calculation of model fitting metrics like loo or waic)

#' @examples 
#' ## Example 1: Fit beta regression models with different link functions to independent data
#' # Load a dataset of independent response observations simulated with a logit link function
#' data(data5)
#' # Note: data were simulated with Coef[1]=1 (intercept), Coef[2]=0.8 (slope), phi=5
#'
#' # Run the betareg function
#' result1 = betareg_stan(des=data5[,1], y=data5[,2], iter=1000)
#'
#' # Fit the model again but this time use a probit link function
#' result2 = betareg_stan(des=data5[,1], y=data5[,2], link="probit", iter=1000)
#'
#' # Compare posterior parameter estimates
#' result1[[1]]
#' result2[[1]]
#'
#' ## Example 2: Fit beta regression models to a simulated dataset in which the data are non-independent
#' \dontrun{
#' # Load a dataset of non-independent response observations simulated with a logit link function
#' data(data7)
#' Note: data were simulated with Coef[1]=1 (intercept), Coef[2]=0.8 (slope), phi=5
#' # Also load the lineage-pair covariance matrix that arose from those simulations
#' data(sim.cov.pairs)
#'
#' # Run the betareg function
#' result1 = betareg_stan(des=data7[,1], y=data7[,2], cov=sim.cov.pairs, iter=1000)
#'
#' # Fit the model again but this time without the covariance matrix
#' result2 = betareg_stan(des=data7[,1], y=data7[,2], iter=1000)
#'
#' # Fit the model a third time with the covariance matrix but now with a probit link function
#' result3 = betareg_stan(des=data7[,1], y=data7[,2], cov=sim.cov.pairs, link="probit", iter=1000)
#'
#' # Compare posterior parameter estimates
#' result1[[1]]
#' result2[[1]]
#' result3[[1]]
#'
#' # Compare the fit of the three models via leave-on-out (loo) cross validation.
#' loo1 = suppressWarnings(loo::loo(result1[[2]]))
#' loo2 = suppressWarnings(loo::loo(result2[[2]]))
#' loo3 = suppressWarnings(loo::loo(result3[[2]]))
#' loo::loo_compare(loo1, loo2, loo3)
#' }

#' @export
betareg_stan = function(des, y, link="logit", covmat=NULL, iter=6000, chains=4, coef.u=0, 
  coef.sd=10, phi.shape=0.01, phi.rate=0.01, physig.u=-1, physig.sd=1, cores=4, ...) {
  # Check response variable range
  if(any(y>1) | any(y<0)) stop("An unbounded response variable has been provided; beta regression not appropriate")
  # Add an intercept column to design matrix if it's not already present
  if(!all(as.matrix(des)[,1]==1)) des=cbind(rep(1, length(y)), des)
  # Adjust response variable to avoid 0 and 1 
  if(0%in%y | 1%in% y) y = y*(length(y) - 1 + 0.5)/length(y)
  # Set link function choice
  link_choice = switch(link, logit = 1, probit = 2, cloglog = 3, loglog = 4, stop("Invalid link function"))
  # Prep stan data
  stan.dat = list(N=length(y), K=ncol(des), Y=y, X=des, coef_mean=coef.u, coef_sd=coef.sd, phi_shape=phi.shape, phi_rate=phi.rate, sig2_mean=physig.u, sig2_sd=physig.sd, link_choice=link_choice)
  # Add covariance to stan data list if provided
  if(!is.null(covmat)) stan.dat$Cp = covmat
  # Choose the model to fit
  stan_model = if(!is.null(covmat)) stanmodels$betamm else stanmodels$betaols
  # Fit the model
  fit = sampling(object = stan_model, data = stan.dat, iter = iter, chains = chains, cores = cores, ...)
  # Get the parameters from the summary
  pars = if(!is.null(covmat)) rstan::summary(fit)$summary[c(1:4,dim(fit)[3]),] else pars = rstan::summary(fit)$summary[c(1:3,dim(fit)[3]),]
  ll = loo::extract_log_lik(fit, parameter_name = "loglik")
  res = list(pars=round(pars,2), ll=ll)
  return(res)
}