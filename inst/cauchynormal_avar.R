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
source("inst/cauchynormal_functions.R")
## 

## Gibbs 
lag <- 100
k <- 100
m <- 5*k
natoms_seq <- c(1, 10, 50, 100)
nrep <- 1e3
results.gibbs <-  foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  run <- sample_unbiasedvar_reservoir(single_kernel, coupled_kernel, rinit, h = h, k = k, m = m, lag = lag, x_0 = x_0, natoms = natoms_seq)
  data.frame(sampler = "gibbs",
             natoms = natoms_seq,
             rep = irep,
             estimator = run$estimator[1,],
             cost = run$cost,
             pih = mean(run$pih),
             varh = run$varh,
             cost_fishyestimation = run$cost_fishyterms)
}
save(results.gibbs, nrep, natoms_seq, k, m, lag, file = "output/cauchynormal.gibbs.avar.RData")
load(file = "output/cauchynormal.gibbs.avar.RData")
head(results.gibbs)


## MRTH
lag <- 75
k <- 75
m <- 5*k
natoms_seq <- c(1, 10, 50, 100)
nrep <- 1e3
results.mrth <-  foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  run <- sample_unbiasedvar_reservoir(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, h = h, k = k, m = m, lag = lag, x_0 = x_0, natoms = natoms_seq)
  data.frame(sampler = "mrth",
             natoms = natoms_seq,
             rep = irep,
             estimator = run$estimator[1,],
             cost = run$cost,
             pih = mean(run$pih),
             varh = run$varh,
             cost_fishyestimation = run$cost_fishyterms)
}
save(results.mrth, nrep, natoms_seq, k, m, lag, file = "output/cauchynormal.mrth.avar.RData")
load(file = "output/cauchynormal.mrth.avar.RData")
head(results.mrth)

# ggplot(results.mrth, aes(x = estimator, fill = factor(natoms))) + geom_histogram() + facet_wrap(~factor(natoms))

get_mean_with_booterror <- function(v, nb = 1000){
  bootresult <- boot::boot(data = v, statistic = function(v,indices) mean(v[indices]), R = nb)
  paste0("[", 
         paste0(prettyNum(as.numeric(quantile(bootresult$t, probs = c(0.025, 0.975))), digits=0, scientific = F, nsmall = 0), collapse = " - "),
         "]")
}
# print variance
get_var_with_booterror <- function(v, nb = 1000){
  bootresult <- boot::boot(data = v, statistic = function(v,indices) var(v[indices]), R = nb)
  paste0("[", 
         paste0(prettyNum(as.numeric(quantile(bootresult$t, probs = c(0.025, 0.975))), digits=2, scientific = T), collapse = " - "),
         "]")
}
# print efficiency
get_eff_with_booterror <- function(v, c, nb = 1000){
  bootresult <- boot::boot(data = cbind(v,c), statistic = function(x,indices) var(x[indices,1])*mean(x[indices,2]), R = nb)
  paste0("[", 
         paste0(prettyNum(as.numeric(quantile(bootresult$t, probs = c(0.025, 0.975))), digits=2, scientific = T), collapse = " - "),
         "]")
}


load(file = "output/cauchynormal.gibbs.avar.RData")
table_gibbs <- results.gibbs %>% ungroup() %>% group_by(natoms) %>%
  summarise(pestimate = get_mean_with_booterror(estimator), pcost = get_mean_with_booterror(cost), pfishycost = get_mean_with_booterror(cost_fishyestimation), pvariance = get_var_with_booterror(estimator), pefficiency = get_eff_with_booterror(estimator,cost))
table_gibbs <- table_gibbs %>% setNames(c("R", "estimate", "total cost", "fishy cost", 
                                        "variance of estimator", "inefficiency"))
print(table_gibbs, width = Inf)
knitr::kable(table_gibbs, digits = 0, row.names = NA, format = 'latex', escape = FALSE) %>%
  cat(., file = 'output/cauchynormal.gibbs.summary.tex')


load(file = "output/cauchynormal.mrth.avar.RData")

table_mrth <- results.mrth %>% ungroup() %>% group_by(natoms) %>%
  summarise(pestimate = get_mean_with_booterror(estimator), pcost = get_mean_with_booterror(cost), pfishycost = get_mean_with_booterror(cost_fishyestimation), pvariance = get_var_with_booterror(estimator), pefficiency = get_eff_with_booterror(estimator,cost))
table_mrth <- table_mrth %>% setNames(c("R", "estimate", "total cost", "fishy cost", 
                                      "variance of estimator", "inefficiency"))
print(table_mrth, width = Inf)
knitr::kable(table_mrth, digits = 0, row.names = NA, format = 'latex', escape = FALSE) %>%
  cat(., file = 'output/cauchynormal.mrth.summary.tex')


