rm(list = ls())
library(unbiasedpoisson)
setmytheme()
library(doParallel)
library(doRNG)
library(dplyr)
# register parallel cores
registerDoParallel(cores = detectCores()-2)
# set RNG seed
set.seed(2)

## import functions for the state space model
source("inst/binomialssm1d_functions.R")
## 

## convergence to stationary
load(file = "output/binomialssm1d.meetings.RData")
# ##
niterations <- 500
ubounds <- sapply(1:niterations, function(t) tv_upper_bound(unlist(meetingtimes), lag, t))
g_tvbounds <- qplot(x = 1:niterations, y = ubounds, geom = "line") +
  ylab("TV distance") + xlab("t")
g_tvbounds <- g_tvbounds + scale_y_log10(breaks = c(0, 0.001, 0.01, 0.1,1)) + geom_rangeframe()
g_tvbounds

## estimation of asymptotic variance using a bad x_0
lag <- 500
k <- lag
m <- 5*k
natoms_seq <- c(1, 10, 20, 50)
nrep <- 5e2

print(x_0)
results.badx0 <-  foreach(irep = 1:nrep) %dorng% {
  # results.test <-  foreach(irep = 1:nrep) %dorng% {
  result <- list()
  result[[1]] <- sample_coupled_chains_and_fish(single_kernel, coupled_kernel, rinit, h = h, k = k, m = m, lag = lag, x_0 = x_0, natoms = max(natoms_seq))
  result[[2]] <- sample_coupled_chains_and_fish(single_kernel, coupled_kernel, rinit, h = h, k = k, m = m, lag = lag, x_0 = x_0, natoms = max(natoms_seq))
  result
  # }
  # sample_unbiasedvar_reservoir(single_kernel, coupled_kernel, rinit, h = h, k = k, m = m, lag = lag, x_0 = x_0, natoms = natoms_seq)
}
save(results.badx0, x_0, nrep, nparticles, natoms_seq, k, m, lag, file = paste0("output/binomialssm1d.uavar.P", nparticles, ".badanchor.RData"))
load(file = paste0("output/binomialssm1d.uavar.P", nparticles, ".badanchor.RData"))
nrep
# results.badx0[[1]]
##
x_0 <- 0.975
print(x_0)
results.goodx0 <-  foreach(irep = 1:nrep) %dorng% {
  result <- list()
  result[[1]] <- sample_coupled_chains_and_fish(single_kernel, coupled_kernel, rinit, h = h, k = k, m = m, lag = lag, x_0 = x_0, natoms = max(natoms_seq))
  result[[2]] <- sample_coupled_chains_and_fish(single_kernel, coupled_kernel, rinit, h = h, k = k, m = m, lag = lag, x_0 = x_0, natoms = max(natoms_seq))
  result
}
save(results.goodx0, x_0, nrep, nparticles, natoms_seq, k, m, lag, file = paste0("output/binomialssm1d.uavar.P", nparticles, ".goodanchor.RData"))
load(file = paste0("output/binomialssm1d.uavar.P", nparticles, ".goodanchor.RData"))

# load(file = "output/binomialssm1d.uavar.badx0.RData")
# results <- results.badx0
# load(file = "output/binomialssm1d.uavar.goodx0.RData")
# results <- results.goodx0
# 
# cost <- sapply(results, function(x) x$cost[length(natoms_seq)])
# mean(cost)
# cost_fishyterms <- sapply(results, function(x) x$cost_fishyterms[length(natoms_seq)])
# mean(cost_fishyterms)
# 
# estimators <- sapply(results, function(x) x$estimator[1,length(natoms_seq)])
# var(estimators)
# bootresults.center <- boot::boot(estimators, function(v, indices) mean(v[indices]), R = 1000)
# hist(bootresults.center$t[,1])
# 
# for (inatoms in seq_along(natoms_seq)){
#   estimators <- sapply(results, function(x) x$estimator[inatoms])
#   bootresults.center <- boot::boot(estimators, function(v, indices) mean(v[indices]), R = 1000)
#   par(mfrow = c(2,1))
#   print(c(mean(estimators)-1.96*sd(estimators)/sqrt(nrep), mean(estimators)+1.96*sd(estimators)/sqrt(nrep)))
#   hist(estimators, xlab = "estimator", nclass = 50, main = paste0("estimator, R = ", natoms_seq[inatoms]))
#   abline(v = c(mean(estimators)-1.96*sd(estimators)/sqrt(nrep), mean(estimators)+1.96*sd(estimators)/sqrt(nrep)), lwd = 2, col = "black")
#   hist(bootresults.center$t[,1], xlab = "mean estimator", nclass = 50, main = paste0("bootstrapped mean estimator, R = ", natoms_seq[inatoms]))
#   abline(v = c(mean(estimators)-1.96*sd(estimators)/sqrt(nrep), mean(estimators)+1.96*sd(estimators)/sqrt(nrep)), lwd = 2, col = "black")
# }


