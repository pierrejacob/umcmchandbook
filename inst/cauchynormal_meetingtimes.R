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
source("inst/cauchynormal_functions.R")
## 

## Gibbs sampler
lag <- 100
nrep <- 1e4
meetingtime_runs <-  foreach(irep = 1:nrep) %dorng% {
  sample_meetingtime(single_kernel, coupled_kernel, rinit, lag = lag)
}
meetingtimes_gibbs <- sapply(meetingtime_runs, function(x) x$meetingtime)
save(meetingtimes_gibbs, nrep, lag, file = "output/cauchynormal.gibbs.meetings.RData")
load(file = "output/cauchynormal.gibbs.meetings.RData")
niterations <- 140
ubounds_gibbs <- sapply(1:niterations, function(t) tv_upper_bound(unlist(meetingtimes_gibbs), lag, t))
g_tvbounds <- qplot(x = 1:niterations, y = ubounds_gibbs, geom = "line") +
  ylab("TV distance") + xlab("t")
g_tvbounds <- g_tvbounds + scale_y_log10(breaks = c(0, 0.001, 0.01, 0.1,1)) + geom_rangeframe()
g_tvbounds

## MRTH
lag <- 100
nrep <- 1e4
meetingtime_runs <-  foreach(irep = 1:nrep) %dorng% {
  sample_meetingtime(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, lag = lag)
}
meetingtimes_mrth <- sapply(meetingtime_runs, function(x) x$meetingtime)
save(meetingtimes_mrth, nrep, lag, file = "output/cauchynormal.mrth.meetings.RData")
load(file = "output/cauchynormal.mrth.meetings.RData")
niterations <- 100
ubounds_mrth <- sapply(1:niterations, function(t) tv_upper_bound(unlist(meetingtimes_mrth), lag, t))
g_tvbounds <- qplot(x = 1:niterations, y = ubounds_mrth, geom = "line") +
  ylab("TV distance") + xlab("t")
g_tvbounds <- g_tvbounds + scale_y_log10(breaks = c(0, 0.001, 0.01, 0.1,1)) + geom_rangeframe()
g_tvbounds

## two TV upper bounds
niterations <- 150
df <- data.frame(sampler = "gibbs", t = 1:niterations, tvub = sapply(1:niterations, function(t) tv_upper_bound(unlist(meetingtimes_gibbs), lag, t)))
df <- rbind(df,
      data.frame(sampler = "mrth", t = 1:niterations, tvub = sapply(1:niterations, function(t) tv_upper_bound(unlist(meetingtimes_mrth), lag, t))))
## plot 
ggplot(df, aes(x = t, y = tvub, colour = sampler)) + geom_line() + scale_y_log10() + 
  ylab("TV distance") + xlab("iteration") + scale_color_manual(values = c("#005BBB", "#FFD500"))
