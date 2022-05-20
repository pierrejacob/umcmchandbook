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

burnin <- 200
nmcmc <- 1e4
thin <- 50
nchains <- 10
results <- foreach(ichain = 1:nchains, .combine = rbind) %dorng% {
  sum_test <- 0
  sum_testsquare <- 0
  sum_fish <- 0
  sum_testtimesfish <- 0
  countereval <- 0
  counterfish <- 0
  
  # history <- rep(0, nmcmc)
  state <- rinit()
  for (imcmc in 1:nmcmc){
    state <- single_kernel(state)
    test_eval <- h(state$position)
    # history[imcmc] <- test_eval 
    if (imcmc > burnin){
      countereval <- countereval + 1
      sum_test <- sum_test + test_eval
      sum_testsquare <- sum_testsquare + test_eval^2
      if ((imcmc %% thin) == 0){
        cat("estimate fish at time", imcmc, "\n")
        counterfish <- counterfish + 1  
        upfrun <- sample_unbiasedfishy(coupled_kernel, h, state, state_x_0)
        sum_fish <- sum_fish + upfrun$estimator
        sum_testtimesfish <- sum_testtimesfish + upfrun$estimator * test_eval
      }
    }
  }
  pih <- sum_test/countereval
  varh <- sum_testsquare/countereval - (pih)^2
  fishyterms <- 2 * (sum_testtimesfish - pih * sum_fish) / counterfish
  estimator <- -varh + fishyterms
  data.frame(chain = ichain, estimator = estimator, varh = varh, pih = pih, fishyterms = fishyterms)
}
save(results, nmcmc, burnin, thin, nchains, file = "output/highdimreg.avarmcmc.RData")
load(file = "output/highdimreg.avarmcmc.RData")
results
mean(results$estimator)

var(results$varh)
var(results$estimator)
var(results$fishyterms)


# matplot(history, type = 'l')


