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
xseq <- seq(from = -30, to = 40, length.out = 1000)
gpi <- qplot(xseq, sapply(xseq, function(z) targetpdf(z)/normalizingconstant), geom = 'line') + xlab(TeX("$\\theta$")) + ylab(TeX('$\\pi(\\theta)$'))
gpi
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
x_0 <- 0
state_x_0_mrth <- rinit_mrth(x_0)

## trajectory

set.seed(2)
timehorizon <- 2e2
history <- rep(0, timehorizon)
state <- rinit_mrth()
for (timeindex in 1:timehorizon){
  state <- single_kernel_mrth(state)
  history[timeindex] <- state$position
}
gtrace = qplot(x = 1:timehorizon, y = history, geom = "blank", xlab = "Iteration", ylab = "State space") + geom_path()
gtrace = gtrace + ylim(-40, 30)
gtrace

ggsave(filename = "mrth.trace.pdf", plot = gtrace, width = 10, height = 5)

#### 
# 
library(ggridges)

set.seed(1)
nchains <- 50000
niterations <- 200
# chains <- matrix(0, nrow = niterations, ncol = nchains)
chains <- foreach(ichain = 1:nchains, .combine = rbind) %dorng% {
  chain <- rep(0, niterations)
  chain_state <- rinit_mrth()
  chain[1] <- chain_state$position
  for (i in 2:niterations){
    chain_state <- single_kernel_mrth(chain_state)
    chain[i] <- chain_state$position
  }
  data.frame(iteration = 1:niterations, c = rep(ichain, niterations), value = chain)
}


names(chains) <- c("iteration", "chain", "value")
chains %>% tail
gridges <- ggplot(chains %>% filter(iteration == 1 | iteration %% 20 == 0), aes(x = value, y = factor(iteration), height = ..density..)) +
  geom_density_ridges(scale = 2, fill = "grey", stat = "density", adjust = .7)
gridges <- gridges + scale_y_discrete(breaks = c(1, 40, 80, 120, 160, 200))
gridges <- gridges + xlab("State space") + ylab("Iteration") + xlim(-40, 30)
gridges <- gridges + coord_flip()
gridges
ggsave(filename = "mrth.ridge.pdf", plot = gridges, width = 10, height = 5)

gtrace = gtrace + scale_x_continuous(breaks = c(1, 40, 80, 120, 160, 200))
gtrace
grt = gridExtra::grid.arrange(gridges, gtrace, nrow = 1)
ggsave(file = "mrth.ridgetrace.pdf", plot = grt, width = 15, height = 5)



