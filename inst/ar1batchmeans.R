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

##### Parallel batch means
# chain = list of chains, each of equal length
# r = 1 is regular BM
# r = 2 is flat-top BM
# r = 3 is lugsail BM
parBM <- function(chains, r = 1){
  # number of chains
  nchains <- length(chains)
  # finding batch size
  b.final <- 0
  for (m in 1:nchains){
    b.final <- b.final + batchSize(chains[[m]])
  }
  b.final <- floor(b.final/nchains)
  n <- length(chains[[1]])
  a <- floor(n/b.final)
  ab <- b.final * a
  trash <- n-ab
  big.chain <- numeric(length = ab*nchains)
  if (ab != n){
    for (i in 1:m){
      big.chain[((i-1)*ab+1):(i*ab)] <- chains[[i]][-(1:trash)]
    }
  } else {
    big.chain <- Reduce("c", chains)
  }
  rtn <- mcse.multi(big.chain, r = r)$cov
  return(rtn)
}
## 
## we will compute the BM estimates with {1,2,4,8} chains and time horizons {1e4, 1e5, 1e6}
timehorizons <- c(1e4, 1e5, 1e6)
nchains <- c(1, 2, 4, 8)
## for each 'nchain' we have 'nrep/nchain' independent copies of the parallel BM estimator
## we also need to discard burn-in, ~500 here
burnin <- 500
load("~/Dropbox/UnbiasedPoissonNumerics/ar1.longruns.RData")
# ##
cat("# independent chains:", nrep, "\n")
cat("time horizon:", timehorizon, "\n")
## add iteration, repetition index
longruns <- longruns %>% group_by(rep) %>% mutate(iteration = row_number()) %>% ungroup()
tail(longruns)
bm.df <- data.frame()
for (itimehorizon in seq_along(timehorizons)){
  print(timehorizons[itimehorizon])
  for (ichain in seq_along(nchains)){
    print(nchains[ichain])
    ngroup <- nrep / nchains[ichain]
    groupstart <- seq(from = 1, to = nrep, by = nchains[ichain])
    bm.df <- rbind(bm.df, foreach (igroup = 1:ngroup, .combine = rbind) %do% {
      bm.df_sub <- data.frame()
      sublongruns <- matrix(0, nrow = timehorizons[itimehorizon]-burnin, ncol = nchains[ichain])
      for (i in 1:nchains[ichain]){
        rep_ <- (groupstart[igroup]:(groupstart[igroup]+nchains[ichain]-1))[i]
        sublongruns[,i] <- longruns[(burnin+1+(rep_-1)*timehorizon):(timehorizons[itimehorizon]+(rep_-1)*timehorizon),] %>% pull(chainh)
      }
      chainlist <- lapply(seq_len(ncol(sublongruns)), function(i) sublongruns[,i])
      for (r in c(1,2,3)){
        bm.df_sub <- rbind(bm.df_sub, data.frame(varestimate = parBM(chainlist, r = r), r = r, group = igroup,
                                         nchains = nchains[ichain], timehorizon = timehorizons[itimehorizon]))
      }
      rm(chainlist)
      rm(sublongruns)
      bm.df_sub
    })
  }
}
tail(bm.df)
save(bm.df, timehorizons, nchains, file = "output/ar1batchmeans.RData")
# 
load(file = "output/ar1batchmeans.RData")
bm.df$nchains <- factor(bm.df$nchains)
head(bm.df)
 
