#' linreg_stan

#' @description The function \code{linreg_stan} fits one of two linear regression models to a dataset in the Stan Bayesian modeling framework via the 'rstan' package. Users can choose to fit either a basic ols model of or a mixed model in which there are covarying residuals. For the latter, users must supply a covariance matrix.

#' @param des A vector of predictor variable observations OR, in the case of multiple predictors, a matrix in which each column is a vector of observations of a given predictor. 'betareg.stan' adds a column of 1s to make this a design matrix whose first column corresponds to the model intercept (unless such a column already exists).
#' @param y A vector of response variable observations. 
#' @param cov Covariance matrix for model residuals (a Cp matrix if analyzing lineage-pair data or a phylogenetic vcv matrix if analyzing bounded species data).
#' @param itnum Number of iterations to run on each chain; defaults to 6000.
#' @param chains Number of chains to run; defaults to 4.
#' @param cores Number of cores to be used; defaults to 4 (one chain per core).
#' @param coef.u Mean of the Gaussian prior for each preditor variable coefficient; defaults to 0.
#' @param coef.sd SD of the Gaussian prior for each preditor variable coefficient; defaults to 10.
#' @param scale.u Mean of the lognormal prior for the scale of the residual covariance; defaults to -1.
#' @param scale.sd SD of the lognormal prior for the scale of the residual covariance; defaults to 1.

#' @return A list containing two elements: (1) the posterior distribution of beta model parameters, and (2) the log-likelihood of the posteriors for use in downstream analyses (e.g. the calculation of model fitting metrics like loo or waic)

#' @examples 
#' ## Example 1: Fit beta regression models with different link functions to independent data
#' # Note: data were simulated with beta[1]=1 (intercept), beta[2]=0.8 (slope), phi=5
#' # Load a data simulated with a logit link function
#' data(data5)
#' # Run the betareg function
#' result1 = betareg_stan(des=data5[,1], y=data5[,2], itnum=1000)

#' # Observe posterior parameter estimates
#' result1[[1]]

#' # Fit the model again but this time without the covariance matrix
#' result2 = betareg_stan(des=data5[,1], y=data5[,2], link="probit", itnum=1000)

#' # Observe posterior parameter estimates
#' result2[[1]]

#' # Compare the fit of the two models via loo and waic
#' loo1 = loo::loo(result1[[2]])
#' loo2 = loo::loo(result2[[2]])
#' waic1 = loo::waic(result1[[2]])
#' waic2 = loo::waic(result2[[2]])
#' loo1
#' loo2
#' waic1
#' waic2
#' loo::loo_compare(loo1, loo2)
#' loo::loo_compare(waic1, waic2)

#' ## Example 2: Fit beta regression models to a simulated dataset in which the data are non-independent
#' # Note: data were simulated with beta[1]=1 (intercept), beta[2]=0.8 (slope), phi=5
#' \dontrun{
#' # Also load the lineage-pair covariance matrix that arose from those simulations
#' data(data7)
#' data(sim.cov.pairs)

#' # Run the betareg function
#' result1 = betareg_stan(des=data7[,1], y=data7[,2], cov=sim.cov.pairs, itnum=1000)

#' # Observe posterior parameter estimates
#' result1[[1]]

#' # Fit the model again but this time without the covariance matrix
#' result2 = betareg_stan(des=data7[,1], y=data7[,2], itnum=1000)

#' # Observe posterior parameter estimates
#' result2[[1]]

#' # Fit the model again with the covariance matrix but now with a probit link function
#' result3 = betareg_stan(des=data7[,1], y=data7[,2], cov=sim.cov.pairs, link="probit", itnum=1000)

#' # Observe posterior parameter estimates
#' result3[[1]]

#' # Compare the fit of the three models via loo
#' loo1 = suppressWarnings(loo::loo(result1[[2]]))
#' loo2 = suppressWarnings(loo::loo(result2[[2]]))
#' loo3 = suppressWarnings(loo::loo(result3[[2]]))
#' loo::loo_compare(loo1, loo2, loo3)
#' }

#' @export
linreg_stan = function(des, y, covmat=NULL, itnum=6000, chains=4, coef.u=0, 
  coef.sd=10, scale.u=-1, scale.sd=1, cores=4, ...) {
  # Add an intercept column to design matrix if it's not already present
  if(!all(as.matrix(des)[,1]==1)) des=cbind(rep(1, length(y)), des)
  # Prep stan data
  stan.dat = list(N=length(y), K=ncol(des), Y=y, X=des, coef_mean=coef.u, coef_sd=coef.sd, sig2_mean=scale.u, sig2_sd=scale.sd)
  # Add covariance to stan data matrix if provided
  if(!is.null(covmat)) stan.dat$Cp = covmat
  # Choose the model to fit
  stan_model = if(!is.null(covmat)) stanmodels$lplmm else stanmodels$ols
  # Fit the model
  fit = sampling(object = stan_model, data = stan.dat, iter = itnum, chains = chains, cores = cores, ...)
  # Get the parameters from the summary
  pars = if(!is.null(covmat)) rstan::summary(fit)$summary[c(1:4,dim(fit)[3]),] else pars = rstan::summary(fit)$summary[c(1:3,dim(fit)[3]),]
  ll = loo::extract_log_lik(fit, parameter_name = "loglik")
  res = list(pars=round(pars,2), ll=ll)
  return(res)
}
