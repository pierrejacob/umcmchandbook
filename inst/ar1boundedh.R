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
source("inst/ar1functions.R")
## bounded test function
h <- function(x) sin(2*pi*x/40)
timehorizon <- 1e5
history <- rep(0, timehorizon)
state <- rinit()
for (timeindex in 1:timehorizon){
  state <- single_kernel(state)
  history[timeindex] <- state$position
}
plot(history, h(history))

## next estimate fishy function
xseq <- seq(from = -30, to = 30, length.out = 50)
nrep <- 1e4
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
head(df)
## plot solution of Poisson equation associated with h
ghtilde <- ggplot(df, aes(x = x, y = estimator)) + geom_hline(yintercept = c(0), linetype = 1, alpha = 0.2) +
  geom_smooth(colour = 'black') + xlab('x') + ylab(TeX("$\\tilde{h}(x)$")) 
ghtilde
ggsave(filename = "output/ar1.boundedh.htilde.pdf", plot = ghtilde, width = 8, height = 5)

## second moment
ghtildesquare <- ggplot(df, aes(x = x, y = estimator^2)) + geom_hline(yintercept = c(0), linetype = 1, alpha = 0.2) +
  geom_smooth(colour = 'black') + xlab('x') + ylab(TeX("$E(\\tilde{H}(x)^2)$")) 
ghtildesquare
ggsave(filename = "output/ar1.boundedh.Htilde.m2.pdf", plot = ghtildesquare, width = 8, height = 5)
