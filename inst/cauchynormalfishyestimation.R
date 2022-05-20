rm(list = ls())
library(unbiasedpoisson)
setmytheme()
library(doParallel)
library(doRNG)
library(dplyr)
# register parallel cores
registerDoParallel(cores = detectCores()-2)
# set RNG seed
set.seed(1)

## import functions for the Cauchy-Normal model
source("inst/cauchynormalfunctions.R")
## 

## obtain 'nrep' unbiased estimators for a grid of value of x
nrep <- 1e3
xseq <- seq(from = -30, to = 40, length.out = 100)
df <- data.frame()
for (x in xseq){
  print(x)
  res_mrth <- foreach(irep = 1:nrep) %dorng% {
    ## start two chains, one at x and one at x_0
    state_x <- rinit_mrth(x)
    state_x_0 <- rinit_mrth(x_0)
    sample_unbiasedfishy(coupled_kernel_mrth, h, state_x, state_x_0)
  }
  res_gibbs <- foreach(irep = 1:nrep) %dorng% {
    ## start two chains, one at x and one at x_0
    state_x <- rinit(x)
    state_x_0 <- rinit(x_0)
    sample_unbiasedfishy(coupled_kernel, h, state_x, state_x_0)
  }

  df <- rbind(df, data.frame(x = x, estimator_gibbs = sapply(res_gibbs, function(x) x$estimator),
                             estimator_mrth = sapply(res_mrth, function(x) x$estimator)))
}
tail(df)
save(df, xseq, nrep, file = "output/cauchynormal.fishyfunction.RData")
load(file = "output/cauchynormal.fishyfunction.RData")
## plot solution of Poisson equation associated with h
ghtilde <- ggplot(df, aes(x = x, y = estimator_gibbs)) + geom_hline(yintercept = c(0), linetype = 1, alpha = 0.2) +
  geom_smooth(span = 10, colour = "#005BBB") + xlab('x') + ylab(TeX("fishy function(x)$")) +
  geom_smooth(aes(y = estimator_mrth), colour = '#FFD500')
ghtilde


## second moment
ghtildesquare <- ggplot(df, aes(x = x, y = estimator_gibbs^2)) + geom_hline(yintercept = c(0), linetype = 1, alpha = 0.2) +
  geom_smooth(span = 10, colour = "#005BBB") + xlab('x') + ylab(TeX("second moment (x)$")) +
  geom_smooth(aes(y = estimator_mrth^2), colour = "#FFD500")
ghtildesquare
