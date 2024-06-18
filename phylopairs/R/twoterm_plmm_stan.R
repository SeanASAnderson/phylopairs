stanrunner.castillo = function(des, y, sp1s, sp2s, mats) {
    stan.dat = list(N=length(y), M=sp1s, P=sp2s, K=ncol(des), Y=y, X=des, C1=mats[[3]], C2=mats[[4]], Z1=mats[[1]], Z2=mats[[2]])
    fit = sampling(
      object = castillo,
      data = stan.dat,
      iter = 6000,
      chains = 4,
      cores=4)
  pars = round(summary(fit)$summary,2)[c(1:4,dim(fit)[3]),]
  ll = loo::extract_log_lik(fit, parameter_name = "loglik")
  res = list(pars=pars, ll=ll)
  return(res)
}
