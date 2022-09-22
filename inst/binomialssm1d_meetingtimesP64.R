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
lag <- 100
# draw meeting times
nrep <- 1e5
nparticles <- 64
meetingtime_runs <- foreach (irep = 1:nrep) %dorng% {
  sample_meetingtime(single_kernel, coupled_kernel, rinit, lag = lag)
}
meetingtimes <- sapply(meetingtime_runs, function(x) x$meetingtime)
hist(meetingtimes-lag)
save(meetingtimes, nrep, lag, nparticles, file = paste0("output/binomialssm1d.meetings.P", nparticles, ".RData"))
load(file = paste0("output/binomialssm1d.meetings.P", nparticles, ".RData"))
##
hist(log(meetingtimes), prob = TRUE)
summary(meetingtimes)

niterations <- 1000
ubounds <- sapply(1:niterations, function(t) tv_upper_bound(unlist(meetingtimes), lag, t))
g_tvbounds <- qplot(x = 1:niterations, y = ubounds, geom = "line") +
  ylab("TV distance") + xlab("t")
g_tvbounds <- g_tvbounds + scale_y_log10(breaks = c(0, 0.001, 0.01, 0.1,1)) + geom_rangeframe()
g_tvbounds
