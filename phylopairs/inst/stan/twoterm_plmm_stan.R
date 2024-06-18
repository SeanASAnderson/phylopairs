#' twoterm_plmm_stan
#'
#' @description Fits the two-term plmm model of Castillo (2007) in a Bayesian
#'   framework. Bayesian sampling is conducted in the Stan software via the 'rstan'  
#'   package. Users supply vectors of observations and an ultrametric phylogenetic tree.
#'   Users can alter parameters for model-parameter prior distributions and MCMC sampling
#'   settings. See details.
#'
#' @details 
#'   The model introduced by Castillo (2007) for analyzing lineage-pair data is an
#'   extension of a phylogenetic linear mixed model (plmm) in which there are two 
#'   random-effect terms, one for the 'species 1' and another for the 'species 2' in 
#'   every pair. The model is \code{Y = Xb + Z1u1 + Z1u2 + eps}, where u1 and u2 are 
#'   vectors of species-specific random  effects that covary according to a phylogenetic 
#'   covariance matrix C. In this model, u1 ~ N(0,physig2*C), u2 ~ N(0,physig2*C), and
#'   eps~N(0, randsig2*I), where randsig2 is sigma^2 scaling parameter for identity 
#'   matrix \code{I} and physig2 is the sigma^2 scaling parameter for phylogenetic 
#'   covariance matrix C. See \code{twoterm_plmm_mats} for details on calculating 
#'   the Z1 and Z2 matrices. 
#'
#'   # Prior Distributions for Model Parameters #
#'   The underlying stan models assume the following prior distributions for model parameters
#'   1. Regression coefficients: Gaussian prior (users can set prior mean and sd)
#'   2. `physig2`: lognormal prior (users can set prior mean and sd)
#'   3. `randsig2`: cauchy prior (users can set location and scale parameters of prior)
#'
#' @param des A vector of predictor variable observations OR, in the case of multiple predictors, a matrix in which each column is a vector of observations of a given predictor. 'betareg.stan' adds a column of 1s to make this a design matrix whose first column corresponds to the model intercept (unless such a column already exists).
#' @param y A vector of response variable observations. 
#' @param model Choice of "ols", "pgls", or "pgls.mm"; defaults to ols.
#' @param cov Covariance matrix for model residuals (a lineage-pair covariance matrix if analyzing lineage-pair data or a phylogenetic vcv matrix if analyzing species data).
#' @param iter Number of iterations to run on each chain; defaults to 2000 (more are often necessary).
#' @param chains Number of MCMC chains to run; defaults to 4.
#' @param cores Number of cores to use. 
#' @param ... additional arguments passed to \code{rstan::sampling}, including control parmeters (see rstan::sampling documentation)
#' @param coef.u Mean of the Gaussian prior for each preditor variable coefficient; defaults to 0.
#' @param coef.sd SD of the Gaussian prior for each preditor variable coefficient; defaults to 10.
#' @param physig2.u Mean of the prior distribution (lognormal) for the scale of the phylogenetic component of residual covariance; defaults to -1.
#' @param physig2.sd SD of the prior distribution (lognormal) for the scale of the phylogenetic component of residual covariance; defaults to 1.
#' @param randsig2.loc Location parameter of the prior distribution (cauchy) for the scale of the independent component of residual covariance ; defaults to 0.
#' @param randsig2.sc Scale parameter of the prior distribution (cauchy) for the scale of the independent component of residual covariance ; defaults to 2.5.
#'
#' @return A list containing two elements: (1) the posterior distribution of beta model parameters, and (2) the log-likelihood of the posteriors for use in downstream analyses (e.g. the calculation of model fitting metrics like loo or waic)
#'
#' @return A list of four matrices: Z1, Z2, cov1, and cov2, where the latter two describe the covariance among the random effects in u1 and u2.
#'
#' @examples 
#' # Simulate a tree
#' lin.tree = phytools::pbtree(n=20)
#' # Generate lineage pairs as the pairwise combinations of species in the tree
#' lin.pairs = data.frame(t(combn(lin.tree$tip.label,2))); colnames(lin.pairs)=c("sp1", "sp2")
#' # Calculate the matrices
#' mats = twoterm_plmm_mats(sp.pairs=lin.pairs, tree=lin.tree)
#' # Check structure of design matrices
#' sapply(mats, dim)
#' head(mats$Z1[,1:5])
#' head(mats$Z2[,1:5])
#' head(mats$cov2[,1:5])
#' head(mats$cov2[,1:5])
#' # Ensure covariance matrices are valid covariance matrices
#' sapply(mats[3:4], covmat.check)

#' @export
twoterm_plmm_stan = function(des, y, sp1s, sp2s, tree, iter=6000, chains=4, coef.u=0, 
  coef.sd=10, physig2.u=-1, physig2.sd=1, randsig2.loc=0, randsig2.sc=2.5, cores=4, ...) {
    mats = twoterm_plmm_mats.R(sp.pairs=cbind(sp1s, sp2s), tree=tree))
    stan.dat = list(N=length(y), M=sp1s, P=sp2s, K=ncol(des), Y=y, X=des, C1=mats[[3]], C2=mats[[4]], Z1=mats[[1]], Z2=mats[[2]], 
      coef_mean=coef.u, coef_sd=coef.sd, sig2_mean=physig2.u, sig2_sd=physig2.sd, sig2_loc=randsig2.loc, sig2_sc=randsig2.sc)
    fit = sampling(object = stanmodels$ols, data = stan.dat, iter = iter, chains = chaines, cores=cores, ...)
    pars = round(summary(fit)$summary,2)[c(1:4,dim(fit)[3]),]
    ll = loo::extract_log_lik(fit, parameter_name = "loglik")
    res = list(pars=pars, ll=ll)
    return(res)
}
