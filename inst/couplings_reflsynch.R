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

## coupling of Normal(mu1, Sigma) and Normal(mu2, Sigma)

Sigma <- matrix(c(1,0.6, 0.6,1), nrow = 2)
Sigma

mu1 <- c(0,0)
mu2 <- c(2,-2)

xmin <- -2.5
xmax <- 4.5
ymin <- -4.5
ymax <- 2.5


data.grid <- expand.grid(s.1 = seq(xmin, xmax, length.out=200), s.2 = seq(ymin, ymax, length.out=200))
q.samp1 <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = mu1, sigma = Sigma))
q.samp2 <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = mu2, sigma = Sigma))
gcont <- ggplot(q.samp1, aes(x=s.1, y=s.2, z=prob)) + 
  geom_contour(bins = 5, color = 'grey80') +
  geom_contour(data=q.samp2, bins = 5, color = 'grey80') +  xlab(TeX("$x_1$"))+ ylab(TeX("$x_2$")) +
  coord_fixed(xlim = c(xmin, xmax), ylim = c(ymin, ymax), ratio = 1) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
gcont

set.seed(1)
nrep <- 20
Sigma_chol <- chol(Sigma)
Sigma_chol_inv <- solve(Sigma_chol)
refl_df <- data.frame()
for (irep in 1:nrep){
  res_ <- unbiasedpoisson::rmvnorm_reflectionmax(mu1, mu2, Sigma_chol, Sigma_chol_inv)
  x1 <- res_$xy[,1]
  x2 <- res_$xy[,2]
  refl_df <- rbind(refl_df, data.frame(rep = irep, x1_x = x1[1], x1_y = x1[2], x2_x = x2[1], x2_y = x2[2]))
}

## hyper plane separating mus
Delta <- Sigma_chol_inv %*% (mu1 - mu2)
e <-  Delta / sqrt(sum(Delta^2))
midpoint <- (mu1+mu2)/2
m <- - (mu2[1]-mu1[1])/(mu2[2]-mu1[2])
# +    geom_segment(data=data.frame(x = mu1[1], y = mu1[2], xend = mu2[1], yend = mu2[2]),
# aes(x = x, y = y, xend = xend, yend = yend, z = NULL))

  
grefl <- gcont + geom_abline(intercept = -m*midpoint[1]+midpoint[2], slope = m, linetype = 3) 
grefl <- grefl + geom_point(data=refl_df, aes(x = x1_x, y = x1_y, z = NULL), size = 5) + 
  geom_point(data=refl_df, aes(x = x2_x, y = x2_y, z = NULL), color = 'grey40', shape = 15, size = 5) + 
  geom_segment(data=refl_df, aes(x = x1_x, xend = x2_x, y = x1_y, yend = x2_y, z = NULL), linetype = 1)
grefl <- grefl + annotate(geom = "label", x = mu1[1], y = mu1[2], label = TeX("$\\mu_1$"), size = 8) +
  annotate(geom = "label", x = mu2[1], y = mu2[2], label = TeX("$\\mu_2$"), size = 8)
grefl
# 
ggsave(filename = "coupling.refl.pdf", plot = grefl, width = 8, height = 6)

## independent coupling
set.seed(1)
nrep <- 20
indep_df <- data.frame()
for (irep in 1:nrep){
  x1 <- mu1 + Sigma_chol %*% rnorm(2)
  x2 <- mu2 + Sigma_chol %*% rnorm(2)
  indep_df <- rbind(indep_df, data.frame(rep = irep, x1_x = x1[1], x1_y = x1[2], x2_x = x2[1], x2_y = x2[2]))
}


gindep <- gcont + geom_point(data=indep_df, aes(x = x1_x, y = x1_y, z = NULL), size = 5) + 
  geom_point(data=indep_df, aes(x = x2_x, y = x2_y, z = NULL), color = 'grey40', shape = 15, size = 5) + 
  geom_segment(data=indep_df, aes(x = x1_x, xend = x2_x, y = x1_y, yend = x2_y, z = NULL), linetype = 1)
gindep <- gindep + annotate(geom = "label", x = mu1[1], y = mu1[2], label = TeX("$\\mu_1$"), size = 8) +
  annotate(geom = "label", x = mu2[1], y = mu2[2], label = TeX("$\\mu_2$"), size = 8)
gindep
# ggsave(filename = "coupling.indep.pdf", plot = gindep, width = 8, height = 6)
# 


## synchronous coupling
set.seed(1)
nrep <- 20
synch_df <- data.frame()
for (irep in 1:nrep){
  z <- rnorm(2)
  x1 <- mu1 + Sigma_chol %*% z
  x2 <- mu2 + Sigma_chol %*% z
  synch_df <- rbind(synch_df, data.frame(rep = irep, x1_x = x1[1], x1_y = x1[2], x2_x = x2[1], x2_y = x2[2]))
}


gsynch <- gcont + geom_point(data=synch_df, aes(x = x1_x, y = x1_y, z = NULL), size = 5) + 
  geom_point(data=synch_df, aes(x = x2_x, y = x2_y, z = NULL), color = 'grey40', shape = 15, size = 5) + 
  geom_segment(data=synch_df, aes(x = x1_x, xend = x2_x, y = x1_y, yend = x2_y, z = NULL), linetype = 1)
gsynch <- gsynch + annotate(geom = "label", x = mu1[1], y = mu1[2], label = TeX("$\\mu_1$"), size = 8) +
  annotate(geom = "label", x = mu2[1], y = mu2[2], label = TeX("$\\mu_2$"), size = 8)
gsynch
# 

ggsave(filename = "coupling.synch.pdf", plot = gsynch, width = 8, height = 6)
# 



