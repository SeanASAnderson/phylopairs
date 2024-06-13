pkgname <- "phylopairs"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "phylopairs-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('phylopairs')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("betareg_stan")
### * betareg_stan

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: betareg_stan
### Title: betareg_stan
### Aliases: betareg_stan

### ** Examples

## Example 1: Fit beta regression models with different link functions to independent data
# Load a data simulated with a logit link function
data(data5)
# Run the betareg function
result1 = betareg_stan(des=data5[,1], y=data5[,2], itnum=1000)
# Observe posterior parameter estimates
result1[[1]]
# Fit the model again but this time without the covariance matrix
result2 = betareg_stan(des=data5[,1], y=data5[,2], link="probit", itnum=1000)
# Observe posterior parameter estimates
result2[[1]]
# Compare the fit of the two models via loo and waic
loo1 = loo::loo(result1[[2]])
loo2 = loo::loo(result2[[2]])
waic1 = loo::waic(result1[[2]])
waic2 = loo::waic(result2[[2]])
loo1
loo2
waic1
waic2
loo::loo_compare(loo1, loo2)
loo::loo_compare(waic1, waic2)
## Example 2: Fit beta regression models to a simulated dataset in which the data are non-independent
## Not run: 
##D # Also load the lineage-pair covariance matrix that arose from those simulations
##D data(data7)
##D data(sim.cov.pairs)
##D # Run the betareg function
##D result1 = betareg_stan(des=data7[,1], y=data7[,2], cov=sim.cov.pairs, itnum=1000)
##D # Observe posterior parameter estimates
##D result1[[1]]
##D # Fit the model again but this time without the covariance matrix
##D result2 = betareg_stan(des=data7[,1], y=data7[,2], itnum=1000)
##D # Observe posterior parameter estimates
##D result2[[1]]
##D # Fit the model again with the covariance matrix but now with a probit link function
##D result3 = betareg_stan(des=data7[,1], y=data7[,2], cov=sim.cov.pairs, link="probit", itnum=1000)
##D # Observe posterior parameter estimates
##D result3[[1]]
##D # Compare the fit of the three models via loo
##D loo1 = suppressWarnings(loo::loo(result1[[2]]))
##D loo2 = suppressWarnings(loo::loo(result2[[2]]))
##D loo3 = suppressWarnings(loo::loo(result3[[2]]))
##D loo_compare(loo1, loo2, loo3)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("betareg_stan", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
