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
customcolor <- rgb(0.29, 0.62, 0.85)

## import functions for the Cauchy-Normal model
source("inst/ar1functions.R")
## 

## next estimate fishy function
load(file = "output/ar1.fishyfunction.RData")
head(df)
## plot solution of Poisson equation associated with h
ghtilde <- ggplot(df, aes(x = x, y = estimator)) + geom_hline(yintercept = c(0), linetype = 1, alpha = 0.2) +
  geom_smooth(colour = "#005BBB") + xlab('x') + ylab(TeX("fishy function(x)$")) 
ghtilde
ggsave(filename = "output/ar1.htilde.pdf", plot = ghtilde, width = 8, height = 5)

## second moment
ghtildesquare <- ggplot(df, aes(x = x, y = estimator^2)) + geom_hline(yintercept = c(0), linetype = 1, alpha = 0.2) +
  geom_smooth(colour = "#005BBB") + xlab('x') + ylab(TeX("second moment (x)$")) 
ghtildesquare
ggsave(filename = "output/ar1.Htilde.m2.pdf", plot = ghtildesquare, width = 8, height = 5)

load("output/ar1.uavar.RData")
head(results)
print(natoms_seq)

##
efficiencysummary <- results %>% group_by(natoms) %>% summarise(meanavar = mean(estimator), varavar = var(estimator), 
                                                                meancost = mean(cost)) %>%
  mutate(inefficiency = format(varavar * meancost, digits = 2))
efficiencysummary
avar_exact <- 1/(1-phi)^2
cat("exact value:", avar_exact, "\n")


table <- results %>% group_by(natoms) %>% summarise(
                                           estimate = mean(estimator),
                                           twostderror = 2 * sqrt(var(estimator)/nrep),
                                           totalcost = mean(cost),
                                           fishycost = mean(cost_fishyestimation),
                                           variance = var(estimator),
                                           variancepretty = prettyNum(var(estimator), digits = 2, scientific=T))
## average cost for proposed estimator with R = 10
proposed_cost <- table$totalcost[2]
## variance of proposed estimator with R = 10
proposed_variance <- table$variance[2]


table$inefficiency <- prettyNum(table$totalcost * table$variance, digits=2, scientific =T)
table <- table %>% select(-variance) %>% setNames(c("R", "estimate  $v^{\\parallel}(P,\\test)$", "$2\\times \\hat{\\sigma}$", "cost", "fishy cost", 
                              "variance of $\\hat{v}(P,\\test)$", "inefficiency"))
print(table)
knitr::kable(table, digits = 0, row.names = NA, format = 'latex', escape = FALSE) %>%
  cat(., file = 'output/ar1.summary.tex')

##### comparison with standard methods

##  load results from long runs of MCMC
load(file = "output/ar1batchmeans.RData")
head(bm.df)
bm.summary.df <- bm.df %>% filter(nchains >= 4) %>% group_by(r, nchains, timehorizon) %>% 
  summarise(mean_estimate = mean(varestimate), 
            cost = mean(nchains) * mean(timehorizon),
            bias = mean(varestimate) - avar_exact,
            mse = mean((varestimate - avar_exact)^2))
bm.summary.df$mse_proposed <- prettyNum(proposed_variance/(floor(bm.summary.df$cost/proposed_cost)), digits = 2, scientific = T)
bm.summary.df$mse <- prettyNum(bm.summary.df$mse, digits = 2, scientific=T)
bm.summary.df$cost <- prettyNum(bm.summary.df$cost, digits = 2, scientific=T)
bm.summary.df$timehorizon <- prettyNum(bm.summary.df$timehorizon, digits = 2, scientific=T)

bm.summary.df <- bm.summary.df %>% setNames(c("r", "\\# chains", "time horizon", "average estimate", "cost",  "bias", "MSE", "MSE proposed"))
print(bm.summary.df)

knitr::kable(bm.summary.df, digits = 0, row.names = NA, format = 'latex', escape = FALSE) %>%
  cat(., file = 'output/ar1.summary.bm.tex')

load(file = "output/ar1spectralvar.RData")
head(spectral.df)
spectral.summary.df <- spectral.df %>% filter(nchains >= 4) %>% group_by(r, nchains, timehorizon) %>% 
  summarise(mean_estimate = mean(varestimate), 
            cost = mean(nchains) * mean(timehorizon),
            bias = mean(varestimate) - avar_exact,
            mse = mean((varestimate - avar_exact)^2))
spectral.summary.df$mse_proposed <- prettyNum(proposed_variance/(floor(spectral.summary.df$cost/proposed_cost)), digits = 2, scientific = T)
spectral.summary.df$mse <- prettyNum(spectral.summary.df$mse, digits = 2, scientific=T)
spectral.summary.df$cost <- prettyNum(spectral.summary.df$cost, digits = 2, scientific=T)
spectral.summary.df$timehorizon <- prettyNum(spectral.summary.df$timehorizon, digits = 2, scientific=T)

spectral.summary.df <- spectral.summary.df %>% setNames(c("r", "\\# chains", "time horizon", "average estimate", "cost",  "bias", "MSE", "MSE proposed"))
print(spectral.summary.df)

knitr::kable(spectral.summary.df, digits = 0, row.names = NA, format = 'latex', escape = FALSE) %>%
  cat(., file = 'output/ar1.summary.spectral.tex')



