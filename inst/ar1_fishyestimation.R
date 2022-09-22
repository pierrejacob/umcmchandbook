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

## estimate fishy function at a grid of values x in (-30,30)
xseq <- seq(from = -30, to = 30, length.out = 50)
nrep <- 1e2
df <- data.frame()
for (x in xseq){
  print(x)
  res_ <- foreach(irep = 1:nrep) %dorng% {
    ## start two chains, one at x and one at x_0
    state_x <- rinit(x)
    state_x_0 <- rinit(x_0)
    sample_unbiasedfishy(coupled_kernel, h, state_x, state_x_0)
  }
  df <- rbind(df, data.frame(x = x, estimator = sapply(res_, function(x) x$estimator)))
}
save(df, xseq, nrep, file = "output/ar1.fishyfunction.RData")

load(file = "output/ar1.fishyfunction.RData")
head(df)
## plot solution of Poisson equation associated with h
ghtilde <- ggplot(df, aes(x = x, y = estimator)) + geom_hline(yintercept = c(0), linetype = 1, alpha = 0.2) +
  geom_smooth(colour = "#005BBB") + xlab('x') + ylab(TeX("fishy function(x)$")) 
ghtilde

## second moment of the estimator of the fishy function as a function of x
ghtildesquare <- ggplot(df, aes(x = x, y = estimator^2)) + geom_hline(yintercept = c(0), linetype = 1, alpha = 0.2) +
  geom_smooth(colour = "#005BBB") + xlab('x') + ylab(TeX("second moment (x)$")) 
ghtildesquare

