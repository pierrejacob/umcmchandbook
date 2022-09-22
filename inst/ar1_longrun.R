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

## path where the (large) files will be saved
savepath <- "~/Dropbox/UnbiasedPoissonNumerics/"

## import functions for the AR(1) model
source("inst/ar1functions.R")
## 

## run independent chains
nrep <- 400
## of length 
timehorizon <- 1e6
longruns <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  chainhistory <- rep(0, timehorizon)
  state <- rinit()
  for (timeindex in 1:timehorizon){
    state <- single_kernel(state)
    chainhistory[timeindex] <- state$position
  }
  data.frame(rep = irep, timehorizon = timehorizon, chainh = chainhistory)
}
save(longruns, timehorizon, nrep, file = paste0(savepath, "ar1.longruns.RData"))


