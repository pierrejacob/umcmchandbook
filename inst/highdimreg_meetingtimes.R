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

source("inst/highdimreg_functions.R")

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
