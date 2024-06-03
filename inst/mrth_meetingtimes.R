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
lag <- 50
nrep <- 50
coupled_runs <-  foreach(irep = 1:nrep) %dorng% {
  sample_coupled_chains(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, lag = lag, m = 500)
}
meetingtimes_mrth <- sapply(coupled_runs, function(x) x$meetingtime)
run <- coupled_runs[[which.max(meetingtimes_mrth)]]

df = data.frame(time = 0:run$iteration, x = run$samples1)
df$y= c(run$samples2, rep(NA, lag))
g <- ggplot(df, aes(x = time, y = x, colour = 'X')) + geom_path() + geom_path(aes(y = y-.1, colour = 'Y')) + 
  geom_vline(xintercept = run$meetingtime) + # geom_vline(xintercept = run$meetingtime - lag, colour = "grey40") +
  scale_color_manual(name = "", values = c("black", "grey40"), breaks = c("X", "Y"), guide = "none") +
  xlab("Iteration") + ylab("State space")
g <- g + annotate(geom = 'label', x = run$meetingtime, y = -8, label = TeX("$\\tau$"), size = 8, parse = T)
# g <- g + annotate(geom = 'label', x = run$meetingtime-lag, y = -8, label = TeX("$\\tau-L$"), size = 8, parse = T, colour = "grey40")
g <- g + annotate(geom = 'label', x = 20, y = 12, label = TeX("$(X_t)$"), size = 8, parse = T)
g <- g + annotate(geom = 'label', x = 20, y = -5, label = TeX("$(Y_t)$"), size = 8, parse = T, colour = "grey40")
g
ggsave(filename = "mrth.lchains.pdf", plot = g, width = 8, height = 6)

df = data.frame(time = 0:run$iteration, x = run$samples1)
df$y= c(rep(NA, lag),run$samples2)
gshift <- ggplot(df, aes(x = time, y = x, colour = 'X')) + geom_path() + geom_path(aes(y = y-.1, colour = 'Y')) + 
  geom_vline(xintercept = run$meetingtime) + 
  scale_color_manual(name = "", values = c("black", "grey40"), breaks = c("X", "Y"), guide = "none") +
  xlab("Iteration") + ylab("State space")
gshift <- gshift + annotate(geom = 'label', x = run$meetingtime, y = -8, label = TeX("$\\tau$"), size = 8, parse = T)
# gshift <- gshift + annotate(geom = 'label', x = run$meetingtime-lag, y = -8, label = TeX("$\\tau-L$"), size = 8, parse = T, colour = "grey40")
gshift <- gshift + annotate(geom = 'label', x = 20, y = 12, label = TeX("$(X_t)$"), size = 8, parse = T)
gshift <- gshift + annotate(geom = 'label', x = 30, y = -5, label = TeX("$(Y_{t-L})$"), size = 8, parse = T, colour = "grey40")
gshift
ggsave(filename = "mrth.lchainsshift.pdf", plot = gshift, width = 8, height = 6)
# 
#### meeting times
lag <- 50
nrep <- 1e4
meetingtime_runs <-  foreach(irep = 1:nrep) %dorng% {
  sample_meetingtime(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, lag = lag)
}
meetingtimes_mrth <- sapply(meetingtime_runs, function(x) x$meetingtime)
ghistmeet = qplot(x = meetingtimes_mrth - lag, geom = "blank") + geom_histogram(aes(y = ..density..)) + xlab("Meeting time - lag") + ylab("Density")
ghistmeet
ggsave(filename = "mrth.histmeet.pdf", plot = ghistmeet, width = 8, height = 6)
gtraj <- gridExtra::grid.arrange(g, gshift, ghistmeet, nrow = 1)
ggsave(filename = "mrth.laggedchains.pdf", plot = gtraj, width = 15, height = 5)


# set.seed(1)
# lags <- c(10, 20, 50, 500)
# meetingtimes <- list()
# nrep <- 1e4
# for (il in 1:length(lags)){
#   meetingtimes[[il]] <-  as.numeric(foreach(irep = 1:nrep, .combine = c) %dorng% {
#     sample_meetingtime(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, lag = lags[il])$meetingtime
#   })
# }
# 
# tvbounds.df <- data.frame()
# niterations <- 500
# for (il in 1:length(lags)){
#   ubounds_mrth <- sapply(1:niterations, function(t) tv_upper_bound(meetingtimes[[il]], lags[il], t))
#   tvbounds.df <- rbind(tvbounds.df, data.frame(iteration = 1:niterations, ubounds = ubounds_mrth, lag = lags[il]))
# }
# 
# g_tvbounds <- ggplot(tvbounds.df, aes(x = iteration, y = ubounds, group = lag, linetype = factor(lag))) + geom_line() +
#   scale_linetype(name = "lag")
# g_tvbounds <- g_tvbounds + xlab("Iteration") + ylab("Distance to stationarity")
# g_tvbounds <- g_tvbounds + theme(legend.key.width = unit(2, 'cm'))
# g_tvbounds
# ggsave(filename = "mrth.tvboundslinear.pdf", plot = g_tvbounds, width = 8, height = 6)
# 
# g_tvboundslog <- g_tvbounds + scale_y_log10()
# ggsave(filename = "mrth.tvboundslog.pdf", plot = g_tvboundslog, width = 8, height = 6)
# 
# g = gridExtra::grid.arrange(g_tvbounds, g_tvboundslog, nrow = 1)
# ggsave(filename = "mrth.tvbounds.pdf", plot = g, width = 15, height = 5)

