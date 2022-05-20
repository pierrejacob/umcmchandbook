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
source("inst/cauchynormalfunctions.R")
## 

targetpdf <- function(z) sapply(z, function(t) exp(unnormalizedlogpdf(t)))
normalizingconstant <- integrate(f = targetpdf,lower = -30, upper = 40)$val

xseq <- seq(from = -30, to = 40, length.out = 100)

gpi <- qplot(xseq, sapply(xseq, function(z) targetpdf(z)/normalizingconstant), geom = 'line') + xlab('x') + ylab(TeX('$\\pi(x)$'))
ggsave(filename = 'output/cauchynormal.target.pdf', plot = gpi, width = 8, height = 4)
gh <- qplot(xseq, sapply(xseq, h), geom = 'line') + xlab('x') + ylab('h(x)')
ggsave(filename = 'output/cauchynormal.h.pdf', plot = gh, width = 8, height = 4)

load(file = "output/cauchynormal.fishyfunction.RData")

## plot solution of Poisson equation associated with h
ghtilde <- ggplot(df, aes(x = x, y = estimator_gibbs)) + geom_hline(yintercept = c(0), linetype = 1, alpha = 0.2) +
  geom_smooth(span = 10, colour = "#005BBB") + xlab('x') + ylab(TeX("fishy function(x)$")) +
  geom_smooth(aes(y = estimator_mrth), colour = '#FFD500')
ghtilde

ggsave(filename = 'output/cauchynormal.htilde.pdf', plot = ghtilde, width = 8, height = 5)

## second moment
ghtildesquare <- ggplot(df, aes(x = x, y = estimator_gibbs^2)) + geom_hline(yintercept = c(0), linetype = 1, alpha = 0.2) +
  geom_smooth(span = 10, colour = "#005BBB") + xlab('x') + ylab(TeX("second moment (x)$")) +
  geom_smooth(aes(y = estimator_mrth^2), colour = "#FFD500")
ghtildesquare
ggsave(filename = 'output/cauchynormal.Htilde.m2.pdf', plot = ghtildesquare, width = 8, height = 5)

head(df)
dfm2 <- df %>% group_by(x) %>% summarise(m2_gibbs = mean(estimator_gibbs^2), m2_mrth = mean(estimator_mrth^2))
gm2 <- ggplot(dfm2, 
       aes(x = x, y = m2_gibbs)) + geom_point(col = 'blue') + geom_point(aes(y = m2_mrth), col = 'red')
gm2
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
xseq_larger <- seq(from = -50, to = 50, length.out = 100)
plot(xseq_larger, predict_m2_gibbs_vectorized(xseq_larger))


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
plot(xseq_larger, predict_m2_mrth_vectorized(xseq_larger))

k <- 100
m <- 500
lag <- 100
natoms <- 5

### now try to use optimal selection probabilities in estimation of asymptotic variance
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
  
  xis <- sqrt( (cc1_df$weight * (sapply(cc1_df$atom.1, h) - pi_h2_))^2 * (predict_m2(cc1_df$atom.1)))
  xis <- xis / sum(xis)
  # xis <- rep(1/length(cc1_df$weight),length(cc1_df$weight))
  # xis <- abs(cc1_df$weight)
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
  # xis <- abs(cc2_df$weight)
  # xis <- rep(1/length(cc2_df$weight),length(cc2_df$weight))
  xis <- sqrt( (cc2_df$weight * (sapply(cc2_df$atom.1, h) - pi_h1_))^2 * (predict_m2(cc2_df$atom.1)))
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

nrep <- 1e3

uavar_opt_runs_gibbs <-  foreach(irep = 1:nrep) %dorng% {
  sample_unbiasedvar_opt(single_kernel, coupled_kernel, rinit, h = h, predict_m2 = predict_m2_gibbs_vectorized, 
                         k = k, m = m, lag = lag, x_0 = x_0, natoms = natoms)
}

uavar_runs_gibbs <-  foreach(irep = 1:nrep) %dorng% {
  sample_unbiasedvar(single_kernel, coupled_kernel, rinit, h = h, 
                         k = k, m = m, lag = lag, x_0 = x_0, natoms = natoms)
}

estimators_gibbs <- sapply(uavar_runs_gibbs, function(x) x$estimator)
cost_gibbs <- sapply(uavar_runs_gibbs, function(x) x$cost)
estimators_opt_gibbs <- sapply(uavar_opt_runs_gibbs, function(x) x$estimator)
cost_opt_gibbs <- sapply(uavar_opt_runs_gibbs, function(x) x$cost)
cat("uniform probabilities estimator:", mean(estimators_gibbs), "+-", 2*sd(estimators_gibbs)/sqrt(length(estimators_gibbs)), "\n")
cat("optimal probabilities estimator:", mean(estimators_opt_gibbs), "+-", 2*sd(estimators_opt_gibbs)/sqrt(length(estimators_opt_gibbs)), "\n")

summary.df <- data.frame(algorithm = c("Gibbs", "Gibbs"), xi = c("uniform", "optimal"), 
                         estimate = c(mean(estimators_gibbs), mean(estimators_opt_gibbs)),
                         twostderror = c(2 * sqrt(var(estimators_gibbs)/nrep), 2 * sqrt(var(estimators_opt_gibbs)/nrep)),
                         cost = c(mean(cost_gibbs), mean(cost_opt_gibbs)),
                         variance = c(var(estimators_gibbs), var(estimators_opt_gibbs)))

## we see a real difference
## next, on MRTH

uavar_opt_runs_mrth <-  foreach(irep = 1:nrep) %dorng% {
  sample_unbiasedvar_opt(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, h = h, predict_m2 = predict_m2_mrth_vectorized, 
                         k = k, m = m, lag = lag, x_0 = x_0, natoms = natoms)
}

uavar_runs_mrth <-  foreach(irep = 1:nrep) %dorng% {
  sample_unbiasedvar(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, h = h, 
                     k = k, m = m, lag = lag, x_0 = x_0, natoms = natoms)
}

estimators_opt_mrth <- sapply(uavar_opt_runs_mrth, function(x) x$estimator)
estimators_mrth <- sapply(uavar_runs_mrth, function(x) x$estimator)
cost_mrth <- sapply(uavar_runs_mrth, function(x) x$cost)
cost_opt_mrth <- sapply(uavar_opt_runs_mrth, function(x) x$cost)

var(estimators_opt_mrth)/var(estimators_mrth)
## the difference is even bigger in that case

cat("uniform probabilities estimator:", mean(estimators_mrth), "+-", 2*sd(estimators_mrth)/sqrt(length(estimators_mrth)), "\n")
cat("optimal probabilities estimator:", mean(estimators_opt_mrth), "+-", 2*sd(estimators_opt_mrth)/sqrt(length(estimators_opt_mrth)), "\n")

summary.df <- rbind(summary.df, data.frame(algorithm = c("MRTH", "MRTH"), xi = c("uniform", "optimal"), 
                         estimate = c(mean(estimators_mrth), mean(estimators_opt_mrth)),
                         twostderror = c(2 * sqrt(var(estimators_mrth)/nrep), 2 * sqrt(var(estimators_opt_mrth)/nrep)),
                         cost = c(mean(cost_mrth), mean(cost_opt_mrth)),
                         variance = c(var(estimators_mrth), var(estimators_opt_mrth))))


summary.df$variance <- prettyNum(summary.df$variance, digits = 2, scientific = T)
summary.df <- summary.df %>% setNames(c("algorithm", "selection probabilities $\\xi$", "estimate  $v^{\\parallel}(P,\\test)$", "$2\\times \\hat{\\sigma}$", "cost",  
                              "variance of $\\hat{v}(P,\\test)$"))
print(summary.df)
knitr::kable(summary.df, digits = 0, row.names = NA, format = 'latex', escape = FALSE) %>%
  cat(., file = 'output/cauchynormal.summary.tex')



