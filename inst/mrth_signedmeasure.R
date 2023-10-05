library(umcmchandbook)
library(tidyverse)
setmytheme()
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
# ggsave(filename = "cauchynormal.target.pdf", plot = gpi, width = 8, height = 5)



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
# x_0 <- 0
# state_x_0_mrth <- rinit_mrth(x_0)
# 
# ## trajectory
# 
# timehorizon <- 1e3
# history <- rep(0, timehorizon)
# state <- rinit_mrth()
# for (timeindex in 1:timehorizon){
#   state <- single_kernel_mrth(state)
#   history[timeindex] <- state$position
# }
# gtrace = qplot(x = 1:timehorizon, y = history, geom = "blank", xlab = "iteration", ylab = "state space") + geom_path(colour = "cornflowerblue")
# gtrace
# 
# # ggsave(filename = "cauchynormal.trace.pdf", plot = gtrace, width = 8, height = 5)



#### meeting times

## MRTH
# lag <- 75
# nrep <- 1e3
# meetingtime_runs <-  foreach(irep = 1:nrep) %dorng% {
#   sample_meetingtime(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, lag = lag)
# }
# meetingtimes_mrth <- sapply(meetingtime_runs, function(x) x$meetingtime)
# save(meetingtimes_mrth, nrep, lag, file = "cauchynormal.mrth.meetings.RData")
# load(file = "cauchynormal.mrth.meetings.RData")
# ghistmeet = qplot(x = meetingtimes_mrth - lag, geom = "blank") + geom_histogram(aes(y = ..density..)) + xlab("meeting time - lag") + ylab("density")
# ghistmeet
# # ggsave(filename = "cauchynormal.histmeet.pdf", plot = ghistmeet, width = 8, height = 5)
# 
# 
# 
# 
# niterations <- 75
# ubounds_mrth <- sapply(1:niterations, function(t) tv_upper_bound(unlist(meetingtimes_mrth), lag, t))
# g_tvbounds <- qplot(x = 1:niterations, y = ubounds_mrth, geom = "line") +
#   ylab("TV distance to stationarity") + xlab("iteration")
# g_tvbounds <- g_tvbounds + scale_y_log10(breaks = c(0, 0.001, 0.01, 0.1,1))
# g_tvbounds
# ggsave(filename = "cauchynormal.tvbounds.pdf", plot = g_tvbounds, width = 8, height = 5)


nrep = 10000
lag = 100
m = 100
cchains <- foreach(irep = 1:nrep) %dorng% {
  umcmchandbook::sample_coupled_chains(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, lag = lag, m = m)
}
names(cchains[[1]])

sapply(cchains, function(x) x$meetingtime)

signedmeasure <- umcmchandbook::c_chains_to_dataframe(cchains, k = 0, m = m, prune = FALSE, dopar = TRUE)
names(signedmeasure) <- c('rep', 'MCMC', 'weight', 'position')
head(signedmeasure)

mcmcmeasure <- signedmeasure %>% filter(MCMC == TRUE) %>% select(-MCMC)
head(mcmcmeasure)

bcmeasure <-  signedmeasure %>% filter(MCMC == FALSE) %>% select(-MCMC)

bcmeasure$position = round(bcmeasure$position, 5)
bcmeasure = bcmeasure %>% arrange(position) %>% select(weight, position) %>% group_by(position) %>%  summarise(weight = sum(weight))
head(bcmeasure)


sum(mcmcmeasure$weight)
sum(bcmeasure %>% filter(weight > 0) %>% pull(weight))
## kernel density estimation
mcmcdensity = density(mcmcmeasure$position, adjust = .5)
biascorrectionpositive = density(x = bcmeasure %>% filter(weight > 0) %>% pull(position), weights = bcmeasure %>% filter(weight > 0) %>% pull(weight), adjust = .5)
biascorrectionnegative = density(x = bcmeasure %>% filter(weight < 0) %>% pull(position), weights = -(bcmeasure %>% filter(weight < 0) %>% pull(weight)), adjust = .5)
## interpolations
fmcmc = approxfun(mcmcdensity$x, mcmcdensity$y, rule = 2)
fbcp = approxfun(biascorrectionpositive$x, biascorrectionpositive$y, rule = 2)
fbcn = approxfun(biascorrectionnegative$x, biascorrectionnegative$y, rule = 2)

curve(fbcp(x), from = -10, to = 20,  ylim = c(-.3, +.3))
curve(-fbcn(x), add=T)
curve(fmcmc(x), add=T, lty = 2)

measureplotdf <- data.frame(x = seq(from = -15, to = +25, length.out = 1000))
measureplotdf$target <- sapply(measureplotdf$x, function(z) targetpdf(z)/normalizingconstant)
measureplotdf$umcmc <- sapply(measureplotdf$x, function(v) fmcmc(v)+fbcp(v)-fbcn(v))
measureplotdf$mcmc <- sapply(measureplotdf$x, function(v) fmcmc(v))
measureplotdf$bc <- sapply(measureplotdf$x, function(v) fbcp(v)-fbcn(v))
measureplotdf$bcp <- sapply(measureplotdf$x, function(v) fbcp(v))
measureplotdf$bcn <- sapply(measureplotdf$x, function(v) fbcn(v))



## unbiased MCMC and target
gsignedtarget = ggplot(measureplotdf, aes(x = x, ymin = 0, ymax = umcmc)) + geom_ribbon(aes(fill = "unbiased measure", color = "unbiased measure")) + 
  geom_line(aes(y = target, color = "target", fill = "target")) + ylim(-1.25*max(measureplotdf$bcn), 1.25*max(measureplotdf$umcmc)) + 
  scale_color_manual(name = "", values = c("black", "grey90"), breaks = c("target", "unbiased measure"), guide = "none") +
  scale_fill_manual(name = "", values = c("black", "grey90"), breaks = c("target", "unbiased measure"))
gsignedtarget <- gsignedtarget + xlab("State space") + ylab("Density")
gsignedtarget
ggsave(filename = "mrth.signedtarget.pdf", plot = gsignedtarget, width = 8, height = 6)


## MCMC and bias correction
gmcmc = ggplot(measureplotdf, aes(x = x, ymin = 0, ymax = mcmc, fill = "MCMC measure")) + geom_ribbon() +  ylim(-1.25*max(measureplotdf$bcn), 1.25*max(measureplotdf$umcmc)) + 
  scale_color_manual(name = "", values = c("black", rgb(0.1,0.1,0.1), "cornflowerblue"), breaks = c("target", "MCMC measure", "bias correction"), guide = "none") +
  scale_fill_manual(name = "", values = c("black", rgb(0.3,0.3,0.3), "cornflowerblue"), breaks = c("target", "MCMC measure", "bias correction")) 
gmcmc <- gmcmc + xlab("State space") + ylab("Density")
gmcmc
ggsave(filename = "mrth.mcmcmeas.pdf", plot = gmcmc, width = 8, height = 6)

gbc = ggplot(measureplotdf, aes(x = x, ymin = 0, ymax = bcp, fill = "positive bias correction")) + geom_ribbon() +  ylim(-1.25*max(measureplotdf$bcn), 1.25*max(measureplotdf$umcmc)) + 
  geom_ribbon(aes(ymin = -bcn, ymax = 0, fill = 'negative bias correction'), alpha = 1) + 
  # geom_line(aes(y = target, fill = "target", color = "target")) +
  scale_color_manual(name = "", values = c("black", "navy", rgb(0.5,0.5,0.5), rgb(0.8,0.8,0.8)), breaks = c("target", "MCMC measure", "positive bias correction", "negative bias correction"), guide = "none") +
  scale_fill_manual(name = "", values = c("black", "navy", rgb(0.5,0.5,0.5), rgb(0.8,0.8,0.8)), breaks = c("target", "MCMC measure", "positive bias correction", "negative bias correction")) 
gbc <- gbc + xlab("State space") + ylab("Density")
gbc

ggsave(filename = "mrth.bcmeas.pdf", plot = gbc, width = 8, height = 6)

gs <- gridExtra::grid.arrange(gsignedtarget, gmcmc, gbc, nrow = 1)
ggsave(filename = "mrth.signedmeasure.pdf", plot = gs, width = 15, height = 5)
