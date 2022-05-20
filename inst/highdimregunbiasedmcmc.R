rm(list = ls())
library(unbiasedpoisson)
setmytheme()
library(doParallel)
library(doRNG)
library(dplyr)
# register parallel cores
registerDoParallel(cores = 10)
# set RNG seed
set.seed(1)

source("inst/highdimregfunctions.R")

lag <- 1e3
k <- lag
m <- 5*k

# result <- sample_unbiasedestimator(single_kernel, coupled_kernel, rinit, h = h, k = k, m = m, lag = lag)
# names(result)
# result$elapsedtime
# result$uestimator

nrep <- 100
results <-  foreach(irep = 1:nrep) %dorng% {
  sample_unbiasedestimator(single_kernel, coupled_kernel, rinit, h = h, k = k, m = m, lag = lag)
}
save(results, nrep, k, m, lag, file = "output/highdimreg.umcmc.RData")
load(file = "output/highdimreg.umcmc.RData")

# plot meeting times and upper bounds on TV
meetingtimes <- sapply(results, function(x) x$meetingtime)
print(summary(meetingtimes-lag))
hist(meetingtimes-lag, prob = TRUE)
niterations <- 1e3
ubounds <- sapply(1:niterations, function(t) tv_upper_bound(unlist(meetingtimes), lag, t))
plot(1:niterations, ubounds, type = "l")


## compute inefficiency
mean_cost <- mean(sapply(results, function(x) x$cost))
mean(sapply(results, function(x) x$uestimator))
variance <- var(sapply(results, function(x) x$uestimator))
cat("inefficiency:", mean_cost * variance, "\n")

## to see if variance is very different because of bias correction term
var(sapply(results, function(x) x$mcmcestimator))
var(sapply(results, function(x) x$correction))
var(sapply(results, function(x) x$uestimator))

