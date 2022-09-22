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

nmcmc <- 1e4
nparticles <- 256
nchains <- detectCores()-2
# #### run long chains
tictoc::tic("bssm 1d long run")
history <- foreach(ichain = 1:nchains, .combine = rbind) %dorng% {
  history_onechain <- matrix(0, nrow = nmcmc, ncol = 1)
  state <- rinit()
  history_onechain[1,] <- state$position
  targetpdf <- rep(0, nmcmc)
  targetpdf[1] <- state$current_target
  for (imcmc in 2:nmcmc){
    state <- single_kernel(state)
    history_onechain[imcmc,] <- state$position
    targetpdf[imcmc] <- state$current_target
  }
  data.frame(chain = ichain, iteration = 1:nmcmc, x1 = history_onechain[,1], targetpdf = targetpdf)
}
elapsed <- tictoc::toc(quiet = T)
tictoc::tic.clear()
print(elapsed$toc-elapsed$tic)
save(nmcmc, nchains, history, file = "output/binomialssm1d.longrun.RData")
load(file = "output/binomialssm1d.longrun.RData")
nmcmc
nchains
# 

ggplot(history %>% filter(iteration >= 1e2), aes(x = iteration, y = x1, group = chain)) + geom_line() + ylab(TeX("$\\alpha$"))
ggplot(history %>% filter(iteration >= 1e2), aes(x = x1)) + geom_histogram() + xlab(TeX("$\\alpha$"))
ggplot(history %>% filter(iteration >= 1e2), aes(x = targetpdf)) + geom_histogram() + xlab("target pdf")
# ##
alltargetpdf <- history %>% pull(targetpdf)
summary(alltargetpdf)
decrorder <- order(unique(alltargetpdf), decreasing = T)
for (i in 1:15){
  print(i)
  hightargetpdf_ <- unique(alltargetpdf)[decrorder[i]]
  print(hightargetpdf_)
  print(dim(history[which(history %>% pull(targetpdf) == hightargetpdf_),])[1])
}
## I bet the chain that reached the highest value got stuck for a while
## let's find that bit of trajectory

history[which(history %>% pull(targetpdf) == unique(alltargetpdf)[decrorder[2]]),]

h <- function(x) x[1]
library(mcmcse) 
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
burnin <- 1000
len <- nmcmc
chainlist <- lapply(1:nchains, function(ichain) history[(nmcmc*(ichain-1)+burnin):(nmcmc*(ichain-1)+len),3])
parBM(chainlist, r = 1)[1,1]
parBM(chainlist, r = 2)[1,1]
parBM(chainlist, r = 3)[1,1]
## BM for each chain
for (ichain in 1:nchains){
  chainlist <- list(history[(nmcmc*(ichain-1)+burnin):(nmcmc*(ichain-1)+len),3])
  print(parBM(chainlist, r = 2)[1,1])
}

