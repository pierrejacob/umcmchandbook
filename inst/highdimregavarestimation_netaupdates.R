rm(list = ls())
library(unbiasedpoisson)
setmytheme()
library(doParallel)
library(doRNG)
library(dplyr)
# register parallel cores
registerDoParallel(cores = 14)
# set RNG seed
set.seed(1)

source("inst/highdimregfunctions.R")
single_kernel <- function(state, nrepeats_eta){
  output <- half_t_kernel(X, X_transpose, y, a0=a0, b0=b0, std_MH=std_MH,
                          state$position$xi, state$position$sigma2,
                          state$position$beta, state$position$eta, 
                          nrepeats_eta = nrepeats_eta,
                          t_dist_df=t_dist_df)
  return(list(position=output))
}

##
coupled_kernel <- function(state1, state2, nrepeats_eta){
  output <- coupled_half_t_kernel(X, X_transpose, y, a0=a0, b0=b0, std_MH=std_MH,
                                  state1$position$xi, state2$position$xi, 
                                  state1$position$sigma2, state2$position$sigma2,
                                  state1$position$beta, state2$position$beta, 
                                  state1$position$eta, state2$position$eta,
                                  nrepeats_eta = nrepeats_eta, t_dist_df=t_dist_df)
  state1 <- list('position'=
                   list('beta'=output$beta_1,
                        'eta'=output$eta_1,
                        'sigma2'=output$sigma2_1,
                        'xi'=output$xi_1))
  state2 <- list('position'=
                   list('beta'=output$beta_2,
                        'eta'=output$eta_2,
                        'sigma2'=output$sigma2_2,
                        'xi'=output$xi_2))
  identical <- ((max(abs(output$beta_1 - output$beta_2)) == 0) &&
                  (max(abs(output$eta_1 - output$eta_2)) == 0)   &&
                  (abs(output$sigma2_1 - output$sigma2_2) == 0)  &&
                  (abs(output$xi_1 - output$xi_2) == 0))
  return(list(state1 = state1, state2 = state2, identical = identical))
}

## convergence to stationary

## draw meeting times
nrep <- 1e3
lag <- 1e3
meetingtime_runs <- data.frame()
nretas <- c(3,5)
for (nreta in nretas){
  print(nreta)
  meetingtime_runs <- rbind(meetingtime_runs, foreach (irep = 1:nrep, .combine = rbind) %dorng% {
    result <- sample_meetingtime(function(s) single_kernel(s, nreta), function(s1, s2) coupled_kernel(s1, s2, nreta), rinit, lag = lag)
    data.frame(nreta = nreta, rep = irep, meeting = result$meetingtime, elapsed = result$elapsed)
  })
  save(meetingtime_runs, nrep, lag, nretas, file = "output/highdimreg.meetings.nreta.RData")
}
##
head(meetingtime_runs)
tail(meetingtime_runs)

load(file = "output/highdimreg.meetings.nreta.RData")
ggplot(meetingtime_runs, aes(x = elapsed, group = factor(nreta))) + geom_histogram() + facet_wrap(~ factor(nreta), ncol = 1)
ggplot(meetingtime_runs, aes(x = meeting, group = factor(nreta))) + geom_histogram() + facet_wrap(~ factor(nreta), ncol = 1)

## compute TV upper bounds
niterations <- 2e3
ubounds.df <- data.frame()
for (nreta_ in nretas){
  meetings <- meetingtime_runs %>% filter(nreta == nreta_) %>% pull(meeting)
  ubounds <- sapply(1:niterations, function(t) tv_upper_bound(meetings, lag, t))
  ubounds.df <- rbind(ubounds.df,
                      data.frame(iteration = 1:niterations,
                                 ubounds = ubounds, 
                                 nreta = nreta_))
}
ggplot(ubounds.df, aes(x = iteration, y = ubounds, linetype = factor(nreta))) + geom_line()

# ## unbiased estimation
lag <- 1e3
k <- lag
m <- 5*k
## R in the paper; number of fishy function estimators per signed measure
natoms <- 5
## number of independent repeats
nrep <- 1e3
#
results <- list()
for (inreta in seq_along(nretas)){
  print(nretas[inreta])
  results[[inreta]] <-  foreach(irep = 1:nrep) %dorng% {
    sample_unbiasedvar_reservoir(function(s) single_kernel(s, nretas[inreta]), 
                                 function(s1, s2) coupled_kernel(s1, s2, nretas[inreta]),
                                 rinit, h = h, k = k, m = m, lag = lag, x_0 = x_0, natoms = natoms)
  }
  save(results, nrep, natoms, k, m, lag, nretas, file = "output/highdimreg.uavar.nreta.RData")
}
# 
load(file = "output/highdimreg.uavar.nreta.RData")
  
length(results)
results[[2]]

names(results[[1]])
results.df <- data.frame()
# nretas <- c(1,3)
for (nreta in nretas){
  results.df <- rbind(results.df, foreach(irep = 1:nrep, .combine=rbind) %do% {
    run <- (results[[nreta]])[[irep]]
    data.frame(nreta = nreta,
               natoms = natoms,
               rep = irep,
               estimator = run$estimator,
               cost = run$cost,
               pih = mean(run$pih),
               varh = run$varh,
               cost_fishyestimation = run$cost_fishyterms,
               fishyterms = run$fishyterms,
               elapsed = run$elapsedtime)
  })
}

head(results.df)
tail(results.df)

# 
table <- results.df %>% group_by(nreta) %>% summarise(
  estimate = mean(estimator),
  twostderror = 2 * sqrt(var(estimator)/nrep),
  elapsed = mean(elapsed),
  totalcost = mean(cost),
  fishycost = mean(cost_fishyestimation),
  variance = var(estimator),
  variancepretty = prettyNum(var(estimator), digits = 2, scientific=T))
table$inefficiency <- prettyNum(table$totalcost * table$variance, digits=2, scientific =T)
print(table)
# ##

### comparison of elapsed times

ggplot(results.df, aes(x = elapsed, group = factor(nreta))) + geom_histogram() + facet_wrap(~nreta, ncol = 1)

ggplot(results.df, aes(x = estimator, group = factor(nreta))) + geom_histogram() + facet_wrap(~nreta, ncol = 1)

ggplot(results.df, aes(y = estimator, x = factor(nreta))) + geom_boxplot() +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "red")

