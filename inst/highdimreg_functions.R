### The code is taken from https://github.com/niloyb/CoupledHalfT
### devtools::install_github("niloyb/CoupledHalfT")
### but is copied in the current package so that the latter contains everything to reproduce
### the figures of the article.
################################################################################
### code to generate synthetic data set
# n <- 100
# p <- 500
# s <- 20
# true_beta <- matrix(0,p,1)
# true_beta[1:s] = 2^(-(seq(s)/4-9/4))
# X <- matrix(rnorm(n*p, mean = 0, sd = 1), nrow = n, ncol = p)
# 
# #Error terms
# error_std <- 2
# error_terms = error_std*rnorm(n, mean = 0, sd = 1)
# y = X%*%true_beta + error_terms
# X_transpose <- t(X)

################################################################################
# # Loading riboflavin dataset
data(riboflavin)

y <- as.vector(riboflavin$y)
X <- as.matrix(riboflavin$x)
colnames(X) <- NULL
rownames(X) <- NULL
X_transpose <- t(X)
n <- length(y)
p <- dim(X)[2]

# Half-t degree of freedom
t_dist_df <- 2
a0 <- 1
b0 <- 1
std_MH <- 0.8

################################################################################
single_kernel <- function(state){
  output <- half_t_kernel(X, X_transpose, y, a0=a0, b0=b0, std_MH=std_MH,
                          state$position$xi, state$position$sigma2,
                          state$position$beta, state$position$eta, 
                          t_dist_df=t_dist_df)
  return(list(position=output))
}

coupled_kernel <- function(state1, state2){
  output <- coupled_half_t_kernel(X, X_transpose, y, a0=a0, b0=b0, std_MH=std_MH,
                                  state1$position$xi, state2$position$xi, state1$position$sigma2,
                                  state2$position$sigma2,
                                  state1$position$beta, state2$position$beta, 
                                  state1$position$eta, state2$position$eta,
                                  t_dist_df=t_dist_df)
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

# initialization from the prior
rinit <- function(position){
  if (missing(position)){
    a0 <- 1
    b0 <- 1
    xi <- (rcauchy(1))^{-2}
    sigma2 <- 1/rgamma(1, shape = a0/2, rate = b0/2)
    eta <- (rt(p, df=t_dist_df))^{-2}
    beta <- rnorm(p)*sqrt(sigma2/(xi*eta))
    return(list(position= list('xi' = xi, 'sigma2' = sigma2, 
                       'beta' = beta, 'eta' = eta)))
  } else {
    return(list(position = position))
  }
}

## test function
h <- function(x){
  return(x$beta[2564])
}
## 2564 is the component that is found by a prelim run (below) to have the largest posterior mean in absolute value
## which here indicates its marginal posterior distribution will be bimodal, so more interesting.

# ## arbitrary point used in the definition of the fishy function
state_x_0 <- rinit()
x_0 <- state_x_0$position
# 
# state_x_0$position$xi
# state_x_0$position$sigma2
# plot(sign(state_x_0$position$beta) * sqrt(abs(state_x_0$position$beta)))
# plot(state_x_0$position$eta, log = 'y')
# range(state_x_0$position$eta)
# 


## unbiased estimation
# lag <- 100
# k <- 100
# m <- 500
# results <- sample_unbiasedvar_reservoir(single_kernel, coupled_kernel, rinit,
#                                         h, k, m, lag, x_0, natoms = 4)
# # names(results)
# results$estimator
# results$fishyterms
# results$cost
# results$cost_fishyterms

# # # MCMC estimate from long run
# mcmciterations <- 5e2
# history <- matrix(0, nrow = mcmciterations, ncol = length(state_x_0$position$beta))
# state <- rinit()
# for (iteration in 1:mcmciterations){
#   state <- single_kernel(state)
#   history[iteration,] <- state$position$beta
# }
# # # 
# posteriormeans <- colMeans(history[100:mcmciterations,])
# plot(posteriormeans)
# abline(v = 2564)
# 
# hist(history[100:mcmciterations,2564])

# hist(history)
# plot(history, type = 'l')
# mean(sapply(unbiased_beta1mean, function(x) x$uestimator))
# sd(sapply(unbiased_beta1mean, function(x) x$uestimator)) /sqrt(nrep)
# mean(history[500:mcmciterations])
# library(doParallel)
# library(doRNG)
# library(dplyr)
# registerDoParallel(cores = 6)
# nrep <- 10
# lag <- 100
# set.seed(1)
# sample_meetingtime(single_kernel, coupled_kernel, rinit, lag = lag)
# 
# meetingtime_runs <- foreach (irep = 1:nrep) %dorng% {
#   sample_meetingtime(single_kernel, coupled_kernel, rinit, lag = lag)
# }
# meetingtimes <- sapply(meetingtime_runs, function(x) x$meetingtime)
# hist(meetingtimes-lag)
# # save(meetingtimes, nrep, lag, file = "output/highdimreg.meetings.RData")
# # load(file = "output/highdimreg.meetings.RData")
# niterations <- 2e2
# ubounds <- sapply(1:niterations, function(t) tv_upper_bound(unlist(meetingtimes), lag, t))
# plot(1:niterations, ubounds)
# 


