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

## convergence to stationary
# draw meeting times
nrep <- 1e3
lag <- 1e3
meetingtime_runs <- foreach (irep = 1:nrep) %dorng% {
  sample_meetingtime(single_kernel, coupled_kernel, rinit, lag = lag)
}
meetingtimes <- sapply(meetingtime_runs, function(x) x$meetingtime)
hist(meetingtimes-lag)
save(meetingtimes, nrep, lag, file = "output/highdimreg.meetings.RData")
load(file = "output/highdimreg.meetings.RData")

niterations <- 2e3
ubounds <- sapply(1:niterations, function(t) tv_upper_bound(unlist(meetingtimes), lag, t))
plot(1:niterations, ubounds, type = "l")

## unbiased estimation
lag <- 1e3
k <- lag
m <- 5*k

natoms_seq <- c(1, 2, 5, 10)

nrep <- 1e3
results <-  foreach(irep = 1:nrep) %dorng% {
  sample_unbiasedvar_reservoir(single_kernel, coupled_kernel, rinit, h = h, k = k, m = m, lag = lag, x_0 = x_0, natoms = natoms_seq)
}
save(results, nrep, natoms_seq, k, m, lag, file = "output/highdimreg.uavar.RData")

load(file = "output/highdimreg.uavar.RData")

names(results[[1]])

results.df <- foreach(irep = 1:nrep, .combine=rbind) %do% {
  run <- results[[irep]]
  data.frame(natoms = natoms_seq,
             rep = irep,
             estimator = run$estimator,
             cost = run$cost,
             pih = mean(run$pih),
             varh = run$varh,
             cost_fishyestimation = run$cost_fishyterms,
             fishyterms = run$fishyterms)
}

table <- results.df %>% group_by(natoms) %>% summarise(
  estimate = mean(estimator),
  twostderror = 2 * sqrt(var(estimator)/nrep),
  totalcost = mean(cost),
  fishycost = mean(cost_fishyestimation),
  variance = var(estimator),
  variancepretty = prettyNum(var(estimator), digits = 2, scientific=T))
table$inefficiency <- prettyNum(table$totalcost * table$variance, digits=2, scientific =T)
print(table)
##



