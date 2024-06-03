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

nsamples <- 10000
# normal means
mu1 <- 1
mu2 <- 2
# std deviation
sigma1 <- 1
sigma2 <- 2


# max coupling
max_samples <- matrix(0, nrow=nsamples, ncol=2)
dp <- function(x) dnorm(x, mean=mu1, sd=sigma1, log=TRUE)
dq <- function(x) dnorm(x, mean=mu2, sd=sigma2, log=TRUE)
rp <- function(n) rnorm(n, mean=mu1, sd=sigma1)
rq <- function(n) rnorm(n, mean=mu2, sd=sigma2)
for (isample in 1:nsamples){
  x <- max_samples[isample, 1] <- rp(1)
  if (dp(x) + log(runif(1)) < dq(x)){
    max_samples[isample, 2] <- x
  } else {
    reject <- TRUE
    y <- NA
    while (reject){
      y <- rq(1)
      reject <- (dq(y) + log(runif(1)) < dp(y))
    }
    max_samples[isample, 2] <- y
  }
}

gscatter <- qplot(x=max_samples[, 1], y=max_samples[, 2], geom="blank") +
  geom_point(alpha=0.15) +
  geom_abline(slope=1, intercept=0, linetype=2)
gscatter <- gscatter + xlab(TeX("$X$")) + ylab("")
gscatter

gmargx <- ggplot(data = data.frame(x=max_samples[, 1])) +
  geom_histogram(aes(x = x, y=..density..)) 
# add the density of the first normal distribution
gmargx <- gmargx + stat_function(fun=function(x) dnorm(x, mu1, sigma1))
# scale_y_reverse()
gmargx <- gmargx
gmargx


gmargy <- ggplot(data = data.frame(x=max_samples[, 2])) +
  geom_histogram(aes(x = x, y=..density..)) 
# add the density of the second normal distribution
gmargy <- gmargy + stat_function(fun=function(x) dnorm(x, mu2, sigma2))
gmargy <- gmargy + scale_y_reverse() + coord_flip()
gmargy

empty <- ggplot()
g <- gridExtra::grid.arrange(empty, gmargx, gmargy, gscatter,
                             ncol=2, nrow=2,
                             widths=c(1, 4), heights=c(1, 4))
g

# load package for ggmarginplot
library(ggExtra)


p <- ggplot(data=data.frame(x = max_samples[,1], y = max_samples[,2]), 
            aes(x = x, y = y)) +
  geom_point() + scale_x_continuous(breaks = c(-1,0,1,2,3)) +
  scale_y_continuous(breaks = seq(from=-2, to = 6, by = 2)) + 
  geom_vline(xintercept = 1, size = 0.5) + geom_hline(yintercept = 2, size = 0.5)

# Plot the scatter plot with marginal histograms
g <- ggMarginal(p, type = "histogram", xparams = list(fill = "grey70"), 
                yparams = list(fill = "grey30"))
g 
ggsave(filename = "coupling.max.png", plot = g, width = 8, height = 6)
ggsave(filename = "coupling.max.pdf", plot = g, width = 8, height = 6)


pp <- function(x) dnorm(x, mean=mu1, sd=sigma1)
qq <- function(x) dnorm(x, mean=mu2, sd=sigma2)

xseq <- seq(from=-5, to=9, length.out=3e3)
pseq <- pp(xseq)
qseq <- qq(xseq)

g1 <- ggplot(data=data.frame(x=xseq, y=pseq), aes(x=x, y=y)) +
  geom_polygon(fill='grey70', alpha=0.85) +
  geom_polygon(aes(y=qseq),
               fill='grey30', alpha=0.85) +
  geom_polygon(aes(y=pmin(pseq,qseq)),
               col='black', alpha=0.5, linetype=2)
g1 <- g1 + xlab("State space") + ylab("Density")
g1 <- g1 + scale_x_continuous(breaks = seq(from=-2, to = 6, by = 1))
# g1 <- g1 + geom_vline(xintercept = c(1,2), size = 0.5)
# g1 <- g1 + annotate(geom = "label", x = mu1, y = -0.03, label = TeX("$\\mu_1$"), fill = "grey70", size = 8)
g1 <- g1 + annotate(geom = "label", x = mu1-1, y = 0.42, label = TeX("$p$"), fill = "grey70", size = 8)
# g1 <- g1 + annotate(geom = "label", x = mu2, y = -0.03, label = TeX("$\\mu_2$"), fill = "grey30", color = "white", size = 8)
g1 <- g1 + annotate(geom = "label", x = mu2+1, y = 0.22, label = TeX("$q$"), fill = "grey30", color = "white", size = 8)
# g1 <- g1 + annotate(geom = "label", x = (mu1+mu2)/2, y = 0.07, label = TeX("$\\min(p,q)$"), size = 8)
g1 <- g1 + ylim(-.03,0.45)
g1



ggsave(filename = "coupling.max.densities.pdf", plot = g1, width = 8, height = 6)


