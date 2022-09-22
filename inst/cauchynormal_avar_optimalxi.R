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
##  use optimal selection probabilities in estimation of asymptotic variance
sample_unbiasedvar_opt <- function(single_kernel, coupled_kernel, rinit, h, predict_m2,
                                   k = 0, m = 1, lag = 1, x_0 = NULL, natoms = 1){
  ## run two coupled chains
  cc1_ <-   unbiasedpoisson::sample_coupled_chains(single_kernel, coupled_kernel, rinit, m = m, lag = lag)
  cc2_ <-   unbiasedpoisson::sample_coupled_chains(single_kernel, coupled_kernel, rinit, m = m, lag = lag)
  cost_umcmc <- cc1_$cost + cc2_$cost
  ## from both signed measures, get unbiased estimator of pi(h) and pi(h^2)
  pi_h1_        <- unbiasedpoisson::H_bar(cc1_, h = h, k = k, m = m)
  pi_h2_        <- unbiasedpoisson::H_bar(cc2_, h = h, k = k, m = m)
  pi_hsquared1_ <- unbiasedpoisson::H_bar(cc1_, h = function(x) h(x)^2, k = k, m = m)
  pi_hsquared2_ <- unbiasedpoisson::H_bar(cc2_, h = function(x) h(x)^2, k = k, m = m)
  ## estimator of pi(h^2) as average of two estimators
  pi_hsquared_ <- (1/2)*(pi_hsquared1_ + pi_hsquared2_)
  ## estimator of pi(h)^2 as product of two independent estimators of pi(h)
  pih_deepbreath_squared_ <- pi_h1_ * pi_h2_
  ## estimate of variance under pi
  varh <- pi_hsquared_ - pih_deepbreath_squared_
  ## estimate of expectation under pi
  pi_h_average <- (1/2)*(pi_h1_ + pi_h2_)
  covariance_term <- 0 
  cost_fishyestimation <- 0
  ## next get fishy function estimator for atoms in cc1_
  ## convert to data frame and prune identical atoms 
  cc1_df <- unbiasedpoisson:::c_chains_to_dataframe(cc1_, k, m, dopar = FALSE, prune = TRUE)
  ## compute optimal selection probabilities
  xis <- sqrt((cc1_df$weight * (sapply(cc1_df$atom.1, h) - pi_h2_))^2 * (predict_m2(cc1_df$atom.1)))
  xis <- xis / sum(xis)
  for (iatom in 1:natoms){
    ## select an atom at random
    index_ell <- sample(x = 1:(dim(cc1_df)[1]), size = 1, prob = xis)
    atom_ell <- cc1_df[index_ell,4:dim(cc1_df)[2]]
    ## start two chains
    state_atom_ <- rinit(as.numeric(atom_ell))
    state_x_0 <- rinit(x_0)
    ## estimate fishy function at that atom
    upfrun <- sample_unbiasedfishy(coupled_kernel, h, state_atom_, state_x_0)
    htilde_atom_ <- upfrun$estimator 
    cost_fishyestimation <- cost_fishyestimation + 2 * upfrun$meetingtime
    ## return asymptotic variance estimator and cost
    covariance_term <- covariance_term + cc1_df$weight[index_ell] * (h(atom_ell) - pi_h2_) * htilde_atom_ / xis[index_ell]
  }
  ## next get fishy function estimator for atoms in cc2_
  ## convert to data frame and prune identical atoms 
  cc2_df <- unbiasedpoisson:::c_chains_to_dataframe(cc2_, k, m, dopar = FALSE, prune = TRUE)
  xis <- sqrt((cc2_df$weight * (sapply(cc2_df$atom.1, h) - pi_h1_))^2 * (predict_m2(cc2_df$atom.1)))
  xis <- xis / sum(xis)
  for (iatom in 1:natoms){
    ## select an atom at random
    index_ell <- sample(x = 1:(dim(cc2_df)[1]), size = 1, prob = xis)
    atom_ell <- cc2_df[index_ell,4:dim(cc2_df)[2]]
    ## start two chains
    state_atom_ <- rinit(as.numeric(atom_ell))
    state_x_0 <- rinit(x_0)
    ## estimate fishy function at that atom
    upfrun <- sample_unbiasedfishy(coupled_kernel, h, state_atom_, state_x_0)
    htilde_atom_ <- upfrun$estimator 
    cost_fishyestimation <- cost_fishyestimation + 2 * upfrun$meetingtime
    ## return asymptotic variance estimator and cost
    covariance_term <- covariance_term + cc2_df$weight[index_ell] * (h(atom_ell) - pi_h1_) * htilde_atom_ / xis[index_ell]
  }
  covariance_term <- covariance_term / natoms
  asymvar_estimator <- covariance_term - varh
  ## returns estimator of the asymptotic variance
  ## estimator of the "covariance" term 
  ## estimator of the variance of h under pi
  ## estimator of pi(h)
  ## and cost measured in number of iterations
  return(list(estimator = asymvar_estimator, covariance_term = covariance_term,
              varh = varh, pih = pi_h_average, cost_fishyestimation = cost_fishyestimation,
              cost_umcmc = cost_umcmc, cost = cost_umcmc + cost_fishyestimation))
}
## function to compute bootstrapped statistics for tables
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


## load fishy function estimates
load(file = "output/cauchynormal.fishyfunction.RData")
head(df)
dfm2 <- df %>% group_by(x) %>% summarise(m2_gibbs = mean(estimator_gibbs^2), m2_mrth = mean(estimator_mrth^2))
gm2 <- ggplot(dfm2, 
              aes(x = x, y = m2_gibbs)) + geom_point(col = 'blue') + geom_point(aes(y = m2_mrth), col = 'red')
gamfit_gibbs <- mgcv::gam(m2_gibbs ~ s(x, bs = "cs"), data = dfm2)
gamfit_mrth <- mgcv::gam(m2_mrth ~ s(x, bs = "cs"), data = dfm2)
gamfit.df <- data.frame(x = xseq)
gamfit.df$predict_gibbs <- as.numeric(mgcv::predict.gam(gamfit_gibbs, gamfit.df, type = "response"))
gamfit.df$predict_mrth <- as.numeric(mgcv::predict.gam(gamfit_mrth, gamfit.df, type = "response"))
gm2 + geom_line(data = gamfit.df, aes(x = x, y = predict_gibbs), col = 'blue') +
  geom_line(data = gamfit.df, aes(x = x, y = predict_mrth),  col = 'red')
predict_m2_gibbs <- function(x){
  if (x <= max(xseq) &&  x>= min(xseq)){
    return(as.numeric(mgcv::predict.gam(gamfit_gibbs, data.frame(x = x), type = "response")))
  } else {
    if (x > max(xseq)){
      return(as.numeric(mgcv::predict.gam(gamfit_gibbs, data.frame(x = max(xseq)), type = "response")))
    } else {
      return(as.numeric(mgcv::predict.gam(gamfit_gibbs, data.frame(x = min(xseq)), type = "response")))
    }
  }
}
predict_m2_gibbs_vectorized <- Vectorize(predict_m2_gibbs)
predict_m2_mrth <- function(x){
  if (x <= max(xseq) &&  x>= min(xseq)){
    return(as.numeric(mgcv::predict.gam(gamfit_mrth, data.frame(x = x), type = "response")))
  } else {
    if (x > max(xseq)){
      return(as.numeric(mgcv::predict.gam(gamfit_mrth, data.frame(x = max(xseq)), type = "response")))
    } else {
      return(as.numeric(mgcv::predict.gam(gamfit_mrth, data.frame(x = min(xseq)), type = "response")))
    }
  }
}
predict_m2_mrth_vectorized <- Vectorize(predict_m2_mrth)

## Gibbs 
lag <- 100
k <- 100
m <- 5*k
natoms <- 10
nrep <- 1e3

results.gibbs.opt <-  foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  run <- sample_unbiasedvar_opt(single_kernel, coupled_kernel, rinit, h = h, predict_m2 = predict_m2_gibbs_vectorized, 
                         k = k, m = m, lag = lag, x_0 = x_0, natoms = natoms)
  data.frame(sampler = "gibbs opt",
             natoms = natoms,
             rep = irep,
             estimator = run$estimator,
             cost = run$cost,
             cost_fishyestimation = run$cost_fishyestimation)
}
save(results.gibbs.opt, nrep, natoms, k, m, lag, file = "output/cauchynormal.gibbs.avar.opt.RData")
load(file = "output/cauchynormal.gibbs.avar.opt.RData")

load(file = "output/cauchynormal.gibbs.avar.RData")
table_gibbs <- results.gibbs %>% filter(natoms == 10) %>% ungroup() %>% 
  summarise(algo = "gibbs",  selection = "uniform", pestimate = get_mean_with_booterror(estimator), pcost = get_mean_with_booterror(cost), pfishycost = get_mean_with_booterror(cost_fishyestimation), pvariance = get_var_with_booterror(estimator), pefficiency = get_eff_with_booterror(estimator,cost))
table_gibbs.opt <- results.gibbs.opt %>% ungroup() %>% 
  summarise(algo = "gibbs", selection = "opt", pestimate = get_mean_with_booterror(estimator), pcost = get_mean_with_booterror(cost), pfishycost = get_mean_with_booterror(cost_fishyestimation), pvariance = get_var_with_booterror(estimator), pefficiency = get_eff_with_booterror(estimator,cost))
rbind(table_gibbs, table_gibbs.opt)

## MRTH
lag <- 75
k <- 75
m <- 5*k
natoms <- 10
results.mrth.opt <-  foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  run <- sample_unbiasedvar_opt(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, h = h, predict_m2 = predict_m2_mrth_vectorized, 
                         k = k, m = m, lag = lag, x_0 = x_0, natoms = natoms)
  data.frame(sampler = "mrth opt",
             natoms = natoms,
             rep = irep,
             estimator = run$estimator,
             cost = run$cost,
             cost_fishyestimation = run$cost_fishyestimation)
  
}
save(results.mrth.opt, nrep, natoms, k, m, lag, file = "output/cauchynormal.mrth.avar.opt.RData")
load(file = "output/cauchynormal.mrth.avar.opt.RData")
table_mrth.opt <- results.mrth.opt %>% ungroup() %>% 
  summarise(algo = "mrth", selection = "opt", pestimate = get_mean_with_booterror(estimator), pcost = get_mean_with_booterror(cost), pfishycost = get_mean_with_booterror(cost_fishyestimation), pvariance = get_var_with_booterror(estimator), pefficiency = get_eff_with_booterror(estimator,cost))
load(file = "output/cauchynormal.mrth.avar.RData")
table_mrth <- results.mrth %>% filter(natoms == 10) %>% ungroup() %>% 
  summarise(algo = "mrth",  selection = "uniform", pestimate = get_mean_with_booterror(estimator), pcost = get_mean_with_booterror(cost), pfishycost = get_mean_with_booterror(cost_fishyestimation), pvariance = get_var_with_booterror(estimator), pefficiency = get_eff_with_booterror(estimator,cost))
rbind(table_mrth, table_mrth.opt)


table_selectionprob <- rbind(table_gibbs, table_gibbs.opt, table_mrth, table_mrth.opt)
table_selectionprob <- table_selectionprob %>% setNames(c("algorithm", "selection probabilities $\\xi$", "estimate", "total cost", "fishy cost", 
                                        "variance of estimator", "inefficiency"))

print(table_selectionprob)
knitr::kable(table_selectionprob, digits = 0, row.names = NA, format = 'latex', escape = FALSE) %>%
  cat(., file = 'output/cauchynormal.selectionprob.tex')


