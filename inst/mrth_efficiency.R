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
x_0 <- 0
state_x_0_mrth <- rinit_mrth(x_0)

#### meeting times
set.seed(1)
nrep <- 1e3
meetingtime_runs <-  foreach(irep = 1:nrep) %dorng% {
  sample_meetingtime(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, lag = 1)
}
meetingtimes_mrth <- sapply(meetingtime_runs, function(x) x$meetingtime)
ghistmeet = qplot(x = meetingtimes_mrth, geom = "blank") + geom_histogram(aes(y = ..density..)) + 
  xlab("meeting time") + ylab("density")
ghistmeet

## test function
h <- function(x){
  return(x)
}


## MRTH asymptotic variance
x_0 <- 0
state_x_0_mrth <- rinit_mrth(x_0)

## unbiased MCMC with given recommendations
L <- floor(quantile(meetingtimes_mrth, probs = .99))
L <- as.numeric(L)
k <- L
m <- 10 * k
natoms_seq <- c(1, 10, 50)
nrep <- 1e3
results.mrth <-  foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  run <- sample_unbiasedvar_reservoir(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, h = h, k = k, m = m, lag = L, x_0 = x_0, natoms = natoms_seq)
  data.frame(sampler = "mrth",
             natoms = natoms_seq,
             rep = irep,
             estimator = run$estimator[1,],
             cost = run$cost,
             pih = mean(run$pih),
             varh = run$varh,
             cost_fishyestimation = run$cost_fishyterms)
}
save(results.mrth, nrep, natoms_seq, k, m, lag, file = "mrth.avar.RData")
load(file = "mrth.avar.RData")
varestimates <- results.mrth %>% filter(natoms == 50) %>% pull(estimator)
ghistavar = qplot(x = varestimates, geom = "blank") + geom_histogram(aes(y = ..density..)) +
  xlab("Asymptotic variance estimates") + ylab("Density")
ghistavar
# ggsave(filename = "mrth.histavar.pdf", plot = ghistavar, width = 8, height = 6)
mean(varestimates) - 1.96 * sd(varestimates)/sqrt(length(varestimates))
mean(varestimates) + 1.96 * sd(varestimates)/sqrt(length(varestimates))
mcmc_asvar <- mean(varestimates)
# costestimates <- results.mrth %>% filter(natoms == 50) %>% pull(cost)
# ghistcostavar = qplot(x = costestimates, geom = "blank") + geom_histogram(aes(y = ..density..)) + xlab("cost of asymptotic variance estimates") + ylab("density")
# ghistcostavar
# ggsave(filename = "mrth.histcostavar.pdf", plot = ghistcostavar, width = 8, height = 5)

##
inef <- function(k, m, L, nrep = 1e3){
  umcmc_runs <-  foreach(irep = 1:nrep) %dorng% {
    sample_unbiasedestimator(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, h = h, k = k, m = m, lag = L)
  }
  umcmc <- sapply(umcmc_runs, function(x) x$uestimator)
  meancost_umcmc <- mean(sapply(umcmc_runs, function(x) x$cost))
  var_umcmc <- var(umcmc)
  return(list(inef = var_umcmc * meancost_umcmc, meancost = meancost_umcmc, var = var_umcmc))
}
##

set.seed(1)
inef.df <- data.frame()
k <- 10
m <- 10*k
L <- 1
result <- inef(k = k, m = m, L = L, nrep = 1e3)
inef.df <- rbind(inef.df, data.frame(k = k, m = m, L = L,
                                     meancost = result$meancost, var = result$var, inef = result$inef))

k <- 10
m <- 10*k
L <- 10
result <- inef(k = k, m = m, L = L, nrep = 1e3)
inef.df <- rbind(inef.df, data.frame(k = k, m = m, L = L,
                                     meancost = result$meancost, var = result$var, inef = result$inef))


k <- 100
m <- 10*k
L <- 1
result <- inef(k = k, m = m, L = L, nrep = 1e3)
inef.df <- rbind(inef.df, data.frame(k = k, m = m, L = L,
                                     meancost = result$meancost, var = result$var, inef = result$inef))


k <- 100
m <- 10*k
L <- 100
result <- inef(k = k, m = m, L = L, nrep = 1e3)
inef.df <- rbind(inef.df, data.frame(k = k, m = m, L = L,
                                     meancost = result$meancost, var = result$var, inef = result$inef))

k <- 400
m <- 10*k
L <- 1
result <- inef(k = k, m = m, L = L, nrep = 1e3)
inef.df <- rbind(inef.df, data.frame(k = k, m = m, L = L,
                                     meancost = result$meancost, var = result$var, inef = result$inef))


k <- 400
m <- 10*k
L <- 400
result <- inef(k = k, m = m, L = L, nrep = 1e3)
inef.df <- rbind(inef.df, data.frame(k = k, m = m, L = L,
                                     meancost = result$meancost, var = result$var, inef = result$inef))

k <- 800
m <- 10*k
L <- 1
result <- inef(k = k, m = m, L = L, nrep = 1e3)
inef.df <- rbind(inef.df, data.frame(k = k, m = m, L = L,
                                     meancost = result$meancost, var = result$var, inef = result$inef))


# ###
save(inef.df, mcmc_asvar, file = "mrth.inef.RData")
load(file = "mrth.inef.RData")
###

ginef <- ggplot(inef.df, aes(x = meancost, y = inef/mcmc_asvar)) + scale_y_log10() +   geom_hline(yintercept = 1) + geom_point(size = 15) +
  geom_label(aes(label = paste0("k=", k, ", ", "L=", L)), parse = FALSE, size = 6) + scale_x_continuous(limits = c(-500, 9e3), breaks = 0:10*1000) + 
  xlab("Mean cost per estimator") + ylab("Inefficiency relative to MCMC") 
ginef

ggsave(filename = "mrth.inef.pdf", plot = ginef, width = 8, height = 6)

