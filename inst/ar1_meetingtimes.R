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
source("inst/ar1_functions.R")
## 
## choose lag
lag <- 500

# draw meeting times
nrep <- 5e3
meetingtime_runs <- foreach (irep = 1:nrep) %dorng% {
  sample_meetingtime(single_kernel, coupled_kernel, rinit, lag = lag)
}
save(meetingtime_runs, nrep, lag, file = "output/ar1.meetingtimes.RData")
load(file = "output/ar1.meetingtimes.RData")

meetingtimes <- sapply(meetingtime_runs, function(x) x$meetingtime)
niterations <- 1000
ubounds <- sapply(1:niterations, function(t) tv_upper_bound(unlist(meetingtimes), lag, t))
g_tvbounds <- qplot(x = 1:niterations, y = ubounds, geom = "line") +
  ylab("TV distance") + xlab("t")
g_tvbounds <- g_tvbounds + scale_y_log10(breaks = c(0, 0.001, 0.01, 0.1,1)) + geom_rangeframe()
g_tvbounds
