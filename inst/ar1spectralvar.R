rm(list = ls())
library(unbiasedpoisson)
setmytheme()
library(doParallel)
# library(doRNG)
library(dplyr)
library(mcmcse) 
# register parallel cores
registerDoParallel(cores = 2)
# set RNG seed
set.seed(1)

source("inst/ar1functions.R")
## 
avar_exact <- 1/(1-phi)^2
library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
sourceCpp("inst/lag.cpp")

##FFT function for calculating RSVe on a centered matrix A
mSVEfft <- function (A, b, method = "bartlett")
{
  n <- nrow(A) # A must be centered matrix
  p <- ncol(A)
  w <- as.vector(lag2(1:b, n = n, b = b, method = method)) # calculate lags
  w <- c(1, w[1:(n-1)], 0, w[(n-1):1])  # starting steps from FFT paper
  w <- Re(fftw_r2c(w))
  FF <- matrix(0, ncol = p, nrow = 2*n)
  FF[1:n,] <- A
  if(p > 1)  # multivariate
  {
    FF <- mvfftw_r2c (FF)
    FF <- FF * matrix(w, nrow = 2*n, ncol = p)
    FF <- mvfftw_c2r(FF) / (2* n )
    return ((t(A) %*% FF[1:n, ]) / n )
  } else if(p == 1)  ##univariate calls
  {
    FF <- fftw_r2c (FF)
    FF <- FF * matrix(w, nrow = 2*n, ncol = p)
    FF <- fftw_c2r(FF) / (2* n )
    return ((t(A) %*% FF[1:n]) / n )
  }
  
}

parSVE <- function(chains, r = 1)
{
  # number of chains
  nchains <- length(chains)
  
  # finding batch size
  b.final <- 0
  for(m in 1:nchains)
  {
    b.final <- b.final + batchSize(chains[[m]])
  }
  b.final <- floor(b.final/nchains)
  
  global.mean <- mean(sapply(chains, mean))
  n <- length(chains)
  
  rsve <- 0
  
  for (m in 1:nchains)
  {
    chain.cen <- scale(chains[[m]], center = global.mean, scale =FALSE)
    foo <- mSVEfft(A = chain.cen, b = b.final, method = "tukey") # change to  "bartlett"
    
    rsve <- rsve + (2*foo - mSVEfft(A = chain.cen, b = floor(b.final/r), method = "tukey"))
  }
  
  rtn <- rsve/nchains
  return(rtn)
}

## 
## we will compute the estimates with {1,2,4,8} chains and time horizons {1e4, 1e5, 1e6}
timehorizons <- c(1e4, 1e5, 1e6)
nchains <- c(1, 2, 4, 8)
## for each 'nchain' we have 'nrep/nchain' independent copies of the parallel BM estimator
## we also need to discard burn-in, ~500 here
burnin <- 500
#### uncomment the following to compute parBM on the long runs
load("~/Dropbox/UnbiasedPoissonNumerics/ar1.longruns.RData")
###
cat("# independent chains:", nrep, "\n")
cat("time horizon:", timehorizon, "\n")
## add iteration, repetition index
longruns <- longruns %>% group_by(rep) %>% mutate(iteration = row_number()) %>% ungroup()
tail(longruns)
spectral.df <- data.frame()
for (itimehorizon in seq_along(timehorizons)){
  print(timehorizons[itimehorizon])
  for (ichain in seq_along(nchains)){
    print(nchains[ichain])
    ngroup <- nrep / nchains[ichain]
    groupstart <- seq(from = 1, to = nrep, by = nchains[ichain])
    spectral.df <- rbind(spectral.df, foreach (igroup = 1:ngroup, .combine = rbind) %do% {
      spectral.df_sub <- data.frame()
      sublongruns <- matrix(0, nrow = timehorizons[itimehorizon]-burnin, ncol = nchains[ichain])
      for (i in 1:nchains[ichain]){
        rep_ <- (groupstart[igroup]:(groupstart[igroup]+nchains[ichain]-1))[i]
        sublongruns[,i] <- longruns[(burnin+1+(rep_-1)*timehorizon):(timehorizons[itimehorizon]+(rep_-1)*timehorizon),] %>% pull(chainh)
      }
      chainlist <- lapply(seq_len(ncol(sublongruns)), function(i) sublongruns[,i])
      for (r in c(1,2,3)){
        spectral.df_sub <- rbind(spectral.df_sub, data.frame(varestimate = parSVE(chainlist, r = r), r = r, group = igroup,
                                         nchains = nchains[ichain], timehorizon = timehorizons[itimehorizon]))
      }
      rm(chainlist)
      rm(sublongruns)
      spectral.df_sub
    })
  }
}
save(spectral.df, timehorizons, nchains, file = "output/ar1spectralvar.RData")
##
load(file = "output/ar1spectralvar.RData")

spectral.df$nchains <- factor(spectral.df$nchains)
spectral.df$r <- factor(spectral.df$r)
tail(spectral.df)


