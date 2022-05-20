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

## import functions for the AR(1) model
source("inst/ar1functions.R")
## 

## choose lag, k and m
k <- 500
m <- 2500
lag <- 250

## number of independent estimators
nrep <- 1e4
## number of atoms at which to estimate fishy function, per signed measure
natoms_seq <- c(1, 10, 50, 100)

# run <- sample_unbiasedvar_reservoir(single_kernel, coupled_kernel, rinit, h = h, k = k, m = m, lag = lag, x_0 = x_0, natoms = natoms_seq)
# run$estimator
# run$cost

## generate independent estimators of the asymptotic variance
results <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  run <- sample_unbiasedvar_reservoir(single_kernel, coupled_kernel, rinit, h = h, k = k, m = m, lag = lag, x_0 = x_0, natoms = natoms_seq)
  data.frame(natoms = natoms_seq,
             rep = irep,
             estimator = run$estimator,
             cost = run$cost,
             pih = mean(run$pih),
             varh = run$varh,
             cost_fishyestimation = run$cost_fishyterms,
             fishyterms = run$fishyterms)
}
save(results, nrep, natoms_seq, k, m, lag, file = "output/ar1.uavar.RData")
load(file = "output/ar1.uavar.RData")
head(results)
# 
results %>% group_by(natoms) %>% summarise(meancost = mean(cost), varestimator = var(estimator)) %>% mutate(inef = meancost * varestimator)

###
results %>% summarise(meanpih = mean(pih), stderror = 2*sd(pih)/sqrt(nrep))
pih_exact <- 0
cat("exact value:", pih_exact, "\n")

###
results %>% summarise(meanvarh = mean(varh), stderror = 2*sd(varh)/sqrt(nrep))
varh_exact <- 1/(1-phi^2)
cat("exact value:", varh_exact, "\n")

###
results %>% group_by(natoms) %>% summarise(meanavar = mean(estimator), varavar = var(estimator), stderror = 2*sd(estimator)/sqrt(nrep), meancost = mean(cost)) %>%
  mutate(inefficiency = format(varavar * meancost, digits = 2)) %>% select(-varavar)

avar_exact <- 1/(1-phi)^2
cat("exact value:", avar_exact, "\n")
