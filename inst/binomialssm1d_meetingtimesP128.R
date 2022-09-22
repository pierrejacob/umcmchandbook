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

## convergence to stationary
lag <- 100
# draw meeting times
nrep <- 1e5
nparticles <- 128
meetingtime_runs <- foreach (irep = 1:nrep) %dorng% {
  sample_meetingtime(single_kernel, coupled_kernel, rinit, lag = lag)
}
meetingtimes <- sapply(meetingtime_runs, function(x) x$meetingtime)
elapsedtimes <- sapply(meetingtime_runs, function(x) x$elapsedtime)
hist(meetingtimes-lag)
hist(elapsedtimes)
save(meetingtimes, elapsedtimes, nrep, lag, nparticles, file = paste0("output/binomialssm1d.meetings.P", nparticles, ".RData"))
load(file = paste0("output/binomialssm1d.meetings.P", nparticles, ".RData"))
##
hist(log(meetingtimes), prob = TRUE)
summary(meetingtimes)

ecdf_f <- ecdf(meetingtimes-lag)
# curve(1-ecdf_f(x), from = 1, to = 5e3, log = 'xy', xlab = 't', ylab = 'P(tau > t)')
g_meetings128 <- ggplot(data=NULL) + stat_function(fun=function(x) 1-ecdf_f(x), n = 1000) + scale_y_log10()+scale_x_log10(limits = c(1, 5e3)) + xlab(TeX("$t$")) + ylab(TeX("$P(\\tau>t)$"))
print(g_meetings128)

x <- seq(from = 200, to = 2000, by = 1)
y <- 1-ecdf_f(x)
lmres <- lm(log(y) ~ log(x))
ggplot(data=data.frame(logx=log(x),logy=log(y)), aes(x = logx, y = logy)) + geom_line() +
  geom_abline(intercept = lmres$coefficients[1], slope = lmres$coefficients[2], colour = 'red') + xlab(TeX("$\\log(t)$")) + ylab(TeX("$\\log P(\\tau>t)$"))
print(lmres)

niterations <- 1000
ubounds <- sapply(1:niterations, function(t) tv_upper_bound(unlist(meetingtimes), lag, t))
g_tvbounds <- qplot(x = 1:niterations, y = ubounds, geom = "line") +
  ylab("TV distance") + xlab("t")
g_tvbounds <- g_tvbounds + scale_y_log10(breaks = c(0, 0.001, 0.01, 0.1,1)) + geom_rangeframe()
g_tvbounds
