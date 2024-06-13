## The following script generates eight sample datasets for illustrative purposes
# Dataset 1: unbounded response, linear rshp between response and predictor, no covariance in residuals
# Dataset 2: unbounded response, no rshp between response and predictor, no covariance in residuals
# Dataset 3: unbounded response, linear rshp between response and predictor, covariance in residuals
# Dataset 4: unbounded response, no rshp between response and predictor, covariance in residuals
# Dataset 5: bounded response, linear rshp between response and predictor on link scale, no covariance in residuals
# Dataset 6: bounded response, no rshp between response and predictor on link scale, no covariance in residuals
# Dataset 7: bounded response, linear rshp between response and predictor on link scale, covariance in residuals
# Dataset 8: bounded response, no rshp between response and predictor on link scale, covariance in residuals

## Load packages
library(phytools)
library(MASS)
library(phylopairs)

## Simulate a tree, find the pairs of taxa, calculate the unscaled lineage-pair covariance matrix
set.seed(4)
# Generate a tree of 20 taxa
tree1=pbtree(n=20)
plot(tree1)
# ensure the tree is ultrametric
is.ultrametric(tree1)
# Find the pairwise combinations of taxa in the tree
taxa.pairs = t(combn(tree1$tip.label,2))
npair = nrow(taxa.pairs)
# Calculate the covariance matrix among the set of pairs from each replicate
cov.pairs = taxapair.vcv(sp.pairs=taxa.pairs, tree=tree1, spnames=TRUE)
# Verify that cov.pairs is a valid covariance matrix
covmat.check(cov.pairs)

### Simulate relationship between predictor and unbounded responses
# Generate a predictor variable
pred=rnorm(nrow(taxa.pairs))
# define an intercept and slope for the relationship between response and predictor
b0=1
b1=0.8
# Generate uncorrelated residuals
sig2.eps=0.5
eps = mvrnorm(mu=rep(0,npair), Sigma=diag(npair)*sig2.eps)
# Generate covarying residuals
eps.cov = mvrnorm(mu=rep(0, npair), Sigma=cov.pairs*sig2.eps)
## Generate Unbounded Responses
lin.nocov = b0 + b1*pred + eps
zero.nocov = b0 + eps
lin.cov = b0 + b1*pred + eps.cov
zero.cov = b0 + eps.cov
## Generate Bounded Responses 
# Note that in this case the relationship is on the scale of the link function
# Assume a logit link function
# Calculate predicted values from inverse of the logit (i.e. the logistic function)
mu.linnocov = exp(lin.nocov)/(1 + exp(lin.nocov)) 
mu.zeronocov = exp(zero.nocov)/(1 + exp(zero.nocov)) 
mu.lincov = exp(lin.cov)/(1 + exp(lin.cov)) 
mu.zerocov = exp(zero.cov)/(1 + exp(zero.cov)) 
# Define precision parameter for beta distribution
phi=5
# convert the mu and phi parameters to alpha/beta parameters used by R to define beta distributed variables
p.par1 = mu.linnocov*phi
p.par2 = mu.zeronocov*phi
p.par3 = mu.lincov*phi
p.par4 = mu.zerocov*phi
q.par1 = (1-mu.linnocov)*phi
q.par2 = (1-mu.zeronocov)*phi
q.par3 = (1-mu.lincov)*phi
q.par4 = (1-mu.zerocov)*phi
beta.linnocov = rbeta(n=npair, shape1 = p.par1, shape2 = q.par1)
beta.zerocov = rbeta(n=npair, shape1 = p.par2, shape2 = q.par2)
beta.lincov = rbeta(n=npair, shape1 = p.par3, shape2 = q.par3)
beta.zeronocov = rbeta(n=npair, shape1 = p.par4, shape2 = q.par4)

# Visualize rshps between response and predictor
# e.g.
plot(pred, lin.cov)
plot(pred, beta.linnocov)

# Join predictor and responses and save datasets
data1 = cbind(pred, lin.nocov)
data2 = cbind(pred, zero.nocov)
data3 = cbind(pred, lin.cov)
data4 = cbind(pred, zero.cov)
data5 = cbind(pred, beta.linnocov)
data6 = cbind(pred, beta.zeronocov)
data7 = cbind(pred, beta.lincov)
data8 = cbind(pred, beta.zerocov)
for(i in 1:8) save(list=paste0("data", i), file=paste0("data/", paste0("data", i), ".rda"))
sim.cov.pairs=cov.pairs
save(sim.cov.pairs, file="sim.cov.pairs")



