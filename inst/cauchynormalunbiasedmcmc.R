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

k <- 100
m <- 2000
lag <- 100
nrep <- 1e3
results <-  foreach(irep = 1:nrep) %dorng% {
  sample_unbiasedestimator(single_kernel, coupled_kernel, rinit, h = h, k = k, m = m, lag = lag)
}
mean_cost <- mean(sapply(results, function(x) x$cost))
mean(sapply(results, function(x) x$uestimator))
variance <- var(sapply(results, function(x) x$uestimator))
cat("inefficiency:", mean_cost * variance, "\n")

results <-  foreach(irep = 1:nrep) %dorng% {
  sample_unbiasedestimator(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, h = h, k = k, m = m, lag = lag)
}
mean_cost <- mean(sapply(results, function(x) x$cost))
mean(sapply(results, function(x) x$uestimator))
variance <- var(sapply(results, function(x) x$uestimator))
cat("inefficiency:", mean_cost * variance, "\n")


natoms <- 10
nrep <- 1e3
results_gibbs <-  foreach(irep = 1:nrep) %dorng% {
  sample_unbiasedvar_reservoir(single_kernel, coupled_kernel, rinit, h = h, k = k, m = m, lag = lag, x_0 = x_0, natoms = natoms)
}

estimators_gibbs <- sapply(results_gibbs, function(x) x$estimator)
mean(estimators_gibbs)
sd(estimators_gibbs)/sqrt(nrep)

results_mrth <-  foreach(irep = 1:nrep) %dorng% {
  sample_unbiasedvar_reservoir(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, h = h, k = k, m = m, lag = lag, x_0 = x_0, natoms = natoms)
}

estimators_mrth <- sapply(results_mrth, function(x) x$estimator)
mean(estimators_mrth)
sd(estimators_mrth)/sqrt(nrep)



