library(umcmchandbook)
library(tidyverse)
setmytheme()
library(latex2exp)
set.seed(1)
library(doParallel)
library(doRNG)
library(dplyr)
# register parallel cores
registerDoParallel(cores = detectCores()-2)

### define kernels
sigma2 <- 100
n <- 3
xobs <- c(-8, 8, 17)

## targetpdf
unnormalizedlogpdf <- function(theta){
  return(-theta^2/(2*sigma2) - sum(log((1+(theta-xobs)^2))))
}

targetpdf <- function(z) sapply(z, function(t) exp(unnormalizedlogpdf(t)))
normalizingconstant <- integrate(f = targetpdf,lower = -30, upper = 40)$val
# xseq <- seq(from = -30, to = 40, length.out = 1000)
# gpi <- qplot(xseq, sapply(xseq, function(z) targetpdf(z)/normalizingconstant), geom = 'line') + xlab(TeX("$\\theta$")) + ylab(TeX('$\\pi(\\theta)$'))
# gpi

## MRTH implementation
target <- function(x){
  if (is.null(dim(x))){
    x <- matrix(x, 1, 1)
  }
  unnormalizedlogpdf(x[1,1])
} 
rinit_mrth <- function(x){
  if (missing(x)){
    x <- rnorm(1, -10, 10)
    return(list(position = matrix(x, 1, 1), current_pdf = target(matrix(x, 1, 1))))
  } else {
    return(list(position = matrix(x, 1, 1), current_pdf = target(matrix(x, 1, 1))))
  }
}
ks <- get_mrth_kernels(target, Sigma_proposal = 2^2)
single_kernel_mrth <- ks$single_kernel
coupled_kernel_mrth <- ks$coupled_kernel


set.seed(1)
nrep <- 1e4

# #### meeting times
# lag <- 50
# meetingtime_runs <-  foreach(irep = 1:nrep) %dorng% {
#   sample_meetingtime(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, lag = lag)
# }
# meetingtimes_mrth <- sapply(meetingtime_runs, function(x) x$meetingtime)
# ghistmeet = qplot(x = meetingtimes_mrth - lag, geom = "blank") + geom_histogram(aes(y = ..density..)) + xlab("Meeting time - lag") + ylab("Density")
# ghistmeet


set.seed(1)
lags <- c(5, 50, 250)
meetingtimes <- list()
for (il in 1:length(lags)){
  meetingtimes[[il]] <-  as.numeric(foreach(irep = 1:nrep, .combine = c) %dorng% {
    sample_meetingtime(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, lag = lags[il])$meetingtime
  })
}

tvbounds.df <- data.frame()
niterations <- 500
for (il in 1:length(lags)){
  ubounds_mrth <- sapply(1:niterations, function(t) tv_upper_bound(meetingtimes[[il]], lags[il], t))
  tvbounds.df <- rbind(tvbounds.df, data.frame(iteration = 1:niterations, ubounds = ubounds_mrth, lag = lags[il]))
}

g_tvbounds <- ggplot(tvbounds.df, aes(x = iteration, y = ubounds, group = lag, linetype = factor(lag))) + geom_line() +
  scale_linetype(name = "lag")
g_tvbounds <- g_tvbounds + xlab("Iteration") + ylab("TV distance to stationarity")
g_tvbounds <- g_tvbounds + theme(legend.key.width = unit(2, 'cm'))
g_tvbounds + geom_hline(yintercept = 1)
# ggsave(filename = "mrth.tvboundslinear.pdf", plot = g_tvbounds, width = 8, height = 6)

g_tvboundslog <- g_tvbounds + geom_hline(yintercept = 1, size = .3) + scale_y_log10() 
g_tvboundslog
ggsave(filename = "mrth.tvboundslog.pdf", plot = g_tvboundslog, width = 8, height = 6)

# g = gridExtra::grid.arrange(g_tvbounds, g_tvboundslog, nrow = 1)
# g
# ggsave(filename = "mrth.tvbounds.pdf", plot = g, width = 15, height = 5)


## 1-Wassserstein bounds
## takes as an input 
## - run: output from sample_coupled_chains
## - tmax: maximum number of iterations at which we want to compute upper bounds
w1bound <- function(run, tmax){
  ubounds <- rep(0, tmax+1)
  # retrieve lag from run
  L <- dim(run$samples1)[1] - dim(run$samples2)[1]
  # iterate from lag up to tau - 1
  for (t in L:(run$meetingtime - 1)){
    ## we are at time t >= lag
    ## find which bounds on |pi_k - pi| to update among k = 0,...,tmax
    ## i.e. find all k such that t = k + j * lag for some j >=  1
    ## i.e. find all k such that k = t - lag - j * lag for some j >= 0
    indices_to_update <- t - (1:ceiling(run$meetingtime / L)) * L
    indices_to_update <- indices_to_update[indices_to_update >= 0 & indices_to_update <= tmax]
    diffchains <- run$samples1[1+t,] - run$samples2[1+t-L,]
    ubounds[1+indices_to_update] <- ubounds[1+indices_to_update] + l2norm(diffchains)
  }
  return(ubounds)
}

niterations <- 500
w1bounds.df <- data.frame()
for (il in 1:length(lags)){
  lag <- lags[il]
  coupled_runs <-  foreach(irep = 1:nrep) %dorng% {
    umcmchandbook::sample_coupled_chains(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, lag = lag, m = 1)
  }
  # ground metric is the Euclidean distance
  l2norm <- function(v) sqrt(sum(v ^ 2))
  w1bounds = foreach(irep = 1:nrep, .combine = '+') %dorng% {
    w1bound(coupled_runs[[irep]], niterations)
  } / nrep
  w1bounds.df <- rbind(w1bounds.df, data.frame(iteration = 0:niterations, ubounds = w1bounds, lag = lag))
}


g_w1bounds <- ggplot(w1bounds.df, aes(x = iteration, y = ubounds, group = lag, linetype = factor(lag))) + geom_line() +
  scale_linetype(name = "lag")
g_w1bounds <- g_w1bounds + xlab("Iteration") + ylab("W1 distance to stationarity")
g_w1bounds <- g_w1bounds + theme(legend.key.width = unit(2, 'cm'))
g_w1boundslog <- g_w1bounds + scale_y_log10() 
g_w1bounds
ggsave(filename = "mrth.w1boundslog.pdf", plot = g_w1boundslog, width = 8, height = 6)




# ## check
# w1bounds_alt <- rep(0, tmax + 1)
# for (time_ in 0:tmax){
#   res_ <- 0
#   for (irep in 1:nrep){
#     run <- coupled_runs[[irep]]
#     L <- dim(run$samples1)[1] - dim(run$samples2)[1]
#     jmax <- floor((run$meetingtime - time_ - 1)/L)
#     if (jmax > 0){
#       for (j in 1:jmax){
#         diffchains <- run$samples1[1 + time_ + j * L,] - run$samples2[1 + time_ + (j-1) * L,]
#         res_ <- res_ + l2norm(diffchains)
#       }
#     }
#   }
#   res_ <- res_ / nrep
#   w1bounds_alt[1 + time_] <- res_
# }
# 
# # head(cbind(w1bounds, w1bounds_alt))
# # summary(abs(w1bounds - w1bounds_alt))
# 
# qplot(x = 0:tmax, y = w1bounds, geom = 'line') + scale_y_log10() + 
#   geom_line(aes(x = 0:tmax, y = w1bounds_alt), color = 'red', linetype = 2)




# ### we can try to obtain a reference by computing the 1-Wasserstein distance between
# ### samples from pi_k and samples from pi_infty where infty is a larger number of iterations (600 below)
# library(transport)
# 
# nchains <- 10000
# niterations <- 600
# # chains <- matrix(0, nrow = niterations, ncol = nchains)
# chains <- foreach(ichain = 1:nchains, .combine = rbind) %dorng% {
#   chain <- rep(0, niterations)
#   chain_state <- rinit_mrth()
#   chain[1] <- chain_state$position
#   for (i in 2:niterations){
#     chain_state <- single_kernel_mrth(chain_state)
#     chain[i] <- chain_state$position
#   }
#   chain
# }
# dim(chains)
# 
# transport::wasserstein1d(chains[,1], chains[,niterations], p = 1)
# 
# subiter <- seq(from = 1, to = niterations, by = 100)
# w1exact <- c()
# for (i in 1:length(subiter)){
#   w1exact <- c(w1exact, transport::wasserstein1d(chains[,subiter[i]], chains[,niterations], p = 1))
# }
# 
# 
# qplot(x = 0:tmax, y = w1bounds, geom = 'line') + scale_y_log10() +
#   geom_line(aes(x = subiter, y = w1exact), color = 'red', linetype = 2)
#   