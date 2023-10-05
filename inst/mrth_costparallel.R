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

## test function
h <- function(x){
  return(x)
}

## unbiased MCMC with given recommendations
L <- 400
k <- L
m <- 10 * k
nrep <- 1e3
umcmc_runs <-  foreach(irep = 1:nrep) %dorng% {
  sample_unbiasedestimator(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, h = h, k = k, m = m, lag = L)
}
save(umcmc_runs, L, k, m, nrep, file = "mrth.costumcmc.RData")
load(file = "mrth.costumcmc.RData")

umcmc_runs[[1]]
umcmc <- sapply(umcmc_runs, function(x) x$uestimator)
cost_umcmc <- sapply(umcmc_runs, function(x) x$cost)
summary(cost_umcmc)
qplot(x = cost_umcmc, geom = "blank") + geom_histogram(aes(y = ..density..)) +
  xlab("cost of each unbiased MCMC estimator") + ylab("density")

nrepeat <- 50
cost_umcmc <- list()
for (irepeat in 1:nrepeat){
  umcmc_runs <-  foreach(irep = 1:nrep) %dorng% {
    sample_unbiasedestimator(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, h = h, k = k, m = m, lag = L)
  }
  cost_umcmc[[irepeat]] <- sapply(umcmc_runs, function(x) x$cost)
}


## 
parallelcomputetime <- function(nproc){
  totaltimes <- rep(0, nrepeat)
  for (irepeat in 1:nrepeat){
    cost_ <- sample(cost_umcmc[[irepeat]], size = length(cost_umcmc[[irepeat]]), replace = FALSE)
    ## wall time of each proc 
    timeproc <- rep(0, nproc)
    ## give cost to idle proc
    for (i in 1:length(cost_)){
      idle = which.min(timeproc)
      timeproc[idle] <- timeproc[idle] + cost_[i]
    }
    totaltimes[irepeat] <- max(timeproc)
  }
  return(totaltimes)
}

save(cost_umcmc, nrepeat, L, k, m, nrep, file = "mrth.costumcmc.RData")
load(file = "mrth.costumcmc.RData")

ctimedf = data.frame()
for (nproc in c(1, 10, 100, 1000, 10000)){
  ctimedf <- rbind(ctimedf, data.frame(np = nproc, ctime = parallelcomputetime(nproc)))
}
gctime <- ggplot(ctimedf %>% group_by(np) %>% summarise(averagectime = mean(ctime)),
       aes(x = np, y = averagectime)) + geom_line() + geom_point() + scale_x_log10() + scale_y_log10(limits = c(1e2, 1e7)) +
  xlab("Number of parallel machines") + ylab("Compute time for 1000 estimates")
gctime

ggsave(filename = "mrth.ctime.pdf", plot = gctime, width = 8, height = 6)

ggplot(ctimedf, aes(x = np, y = ctime, group = np)) + geom_boxplot() + scale_x_log10() + scale_y_log10(limits = c(1, 1e7)) +
  xlab("# parallel machines") + ylab("time to generate 1000 unbiased estimators")

ggplot(ctimedf, aes(x = np, y = ctime, group = np)) + geom_point() + scale_x_log10() + scale_y_log10(limits = c(1e3, 1e7)) +
  xlab("# parallel machines") + ylab("time to generate 1000 unbiased estimators")


