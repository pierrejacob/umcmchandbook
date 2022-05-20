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
source("inst/ar1functions.R")
## 
## choose lag, k and m
k <- 500
m <- 20000
lag <- 500

##
nrep <- 1e4
results <-  foreach(irep = 1:nrep) %dorng% {
  sample_unbiasedestimator(single_kernel, coupled_kernel, rinit, h = h, k = k, m = m, lag = lag)
}
save(results, nrep, k, m, lag, file = "output/ar1.umcmc.RData")
load(file = "output/ar1.umcmc.RData")

# plot meeting times and upper bounds on TV
meetingtimes <- sapply(results, function(x) x$meetingtime)
print(summary(meetingtimes-lag))
hist(meetingtimes-lag, prob = TRUE)
niterations <- 5e2
ubounds <- sapply(1:niterations, function(t) tv_upper_bound(unlist(meetingtimes), lag, t))
plot(1:niterations, ubounds, type = "l")


## compute inefficiency
mean_cost <- mean(sapply(results, function(x) x$cost))
mean(sapply(results, function(x) x$uestimator))
variance <- var(sapply(results, function(x) x$uestimator))
cat("inefficiency:", mean_cost * variance, "\n")
## to be compared with asymptotic variance, which is equal to
1/(1-phi)^2

## to see if variance is very different because of bias correction term
var(sapply(results, function(x) x$mcmcestimator))
var(sapply(results, function(x) x$correction))
var(sapply(results, function(x) x$uestimator))
