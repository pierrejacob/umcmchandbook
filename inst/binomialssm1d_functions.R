library(unbiasedpoisson)
setmytheme()
### load data as in controlled SMC paper
re <- try(read.csv("inst/thaldata.csv", header = FALSE, col.names = FALSE))
if (class(re) == "try-error"){
  print("try again")
  print(getwd())
  re <- try(read.csv("thaldata.csv", header = FALSE, col.names = FALSE))
  head(re)
}
observations <- matrix(re[1:3000,1], ncol = 1)
datalength <- nrow(observations)
### plot data
g <- qplot(x = 1:datalength, y = observations, geom = "line") + xlab("time") + ylab("response")
g

## subset data for faster computation
subobs <- observations[floor(seq.int(from = 1000, to = 1500, length.out = 100)),,drop=F]
g <- qplot(x = 1:length(subobs), y = subobs, geom = "line") + xlab("time") + ylab("response")
g

observations <- subobs
datalength <- length(observations)

### logistic function
logistic <- function(x) 1/(1+exp(-x))

### Original model
## x_0 ~ Normal(0,1)
## x_t | x_{t-1} ~ Normal(alpha x_{t-1}, sigma^2)
## y_t | x_t ~ Binomial(50, logistic(x_t))
### comments
## the latent process is AR(1) with parameters theta = alpha and sigma2 is fixed here
## given the latent process, the logistic function maps to (0,1)
## and the observation y_t are Binomial with 50 trials and probability logistic(x_t)

## we define the twisted model in the terminology of controlled SMC
## all formulas are derived from the controlled SMC paper and its supplementary materials

## the matrix abc is of dimension (T+1) x 3
## the first column is a, second column b, third column c
## at time 0, the initial distribution uses a_0, b_0, c_0 to generate particles (actually, c_0 is not used)
rinit_twisted <- function(nparticles, theta, abc){
  k_0 <- 1/(1+2*abc[1,1])
  return(matrix(rnorm(nparticles, mean = - k_0 * abc[1,2], sd = sqrt(k_0)), ncol = nparticles))
}

sigma2 <- 1.5

## at time t when propagating x_{t-1} to x_t...
## using  a_t, b_t
## the arguments are xparticles  (x_{t-1}), theta (parameter), time (t), abc (contains a_t, b_t)
rtransition_twisted <- function(xparticles, theta, time, abc){
  alpha <- theta[1]
  # sigma2 <- theta[2]
  k_t <- 1/(1/sigma2 + 2*abc[time+1,1])
  mean_ <- k_t * ((alpha / sigma2) * xparticles - abc[time+1,2])
  return(matrix(mean_ + rnorm(nparticles, mean = 0, sd = sqrt(k_t)), nrow = 1))
}

## at time t, when weighting x_t using observation y_t
## the arguments are xparticles  (x_t), theta (parameter), time (t), observation (y_t), abc (contains a_t, b_t, c_t and also a_{t+1}, b_{t+1}, c_{t+1})
potential_twisted <- function(xparticles, theta, time, observation, abc){
  x <- xparticles[1,] # particles x_t
  N <- length(x) # number of particles
  alpha <- theta[1] # explicit names for parameters
  # sigma2 <- theta[2]
  sigma <- sqrt(sigma2)
  a <- abc[time+1,1] # shorted names for a_t, b_t, c_t
  b <- abc[time+1,2]
  c <- abc[time+1,3]
  logpsi_t <- - a * x^2 - b * x - c # log(psi_t(x_t))
  # then there are three cases: t = 0, t >= 1 & t < T, and t = T
  if (time == 0){
    k_0 <- 1/(1+2*a)
    logmu_psi_0 <- 0.5 * log(k_0) + (0.5 * k_0 * (b^2) - c)
    logG_0 <- rep(0, N) # no observation at time zero
    k_1 <- 1/(1/sigma2 + 2*abc[2,1])
    logM1psi1 <- 0.5*log(k_1) - log(sigma) + (0.5 * k_1 * (alpha/sigma2 * x - abc[2,2])^2 - 0.5 * (alpha/sigma*x)^2 - abc[2,3])
    return(logmu_psi_0 + logG_0 + logM1psi1 - logpsi_t)
  } else {
    logG_t <- dbinom(observation, 50, prob = logistic(x), log = TRUE)
    if (time >= 1 && time < datalength){
      k_tp1 <- 1/(1/sigma2 + 2*abc[time+2,1])
      logMtp1psitp1 <- 0.5*log(k_tp1) - log(sigma) + (0.5 * k_tp1 * (alpha/sigma2 * x - abc[time+2,2])^2 - 0.5 * (alpha/sigma*x)^2 - abc[time+2,3])
      return(logG_t + logMtp1psitp1 - logpsi_t)
    } else { # (time == datalength){
      return(logG_t - logpsi_t)
    }
  }
}

particle_filter_twisted <- function(nparticles, theta, observations, abc){
  datalength <- nrow(observations)
  # initialization of particles
  xparticles <- rinit_twisted(nparticles, theta, abc)
  xparticles_all <- array(dim = c(datalength+1, 1, nparticles))
  xparticles_all[1,,] <- xparticles
  # log weights at time zero
  logweights <- potential_twisted(xparticles, theta, 0, NULL, abc)
  logweights_all <- matrix(nrow = datalength+1, ncol = nparticles)
  logweights_all[1,] <- logweights
  maxlogweights <- max(logweights)
  weights <- exp(logweights - maxlogweights)
  # log-likelihood estimator
  ll <- maxlogweights + mean(weights)
  ll_all <- rep(0, datalength+1)
  ll_all[1] <- ll
  # normalize weights
  normweights <- weights / sum(weights)
  # draw ancestors using systematic resampling
  ancestors <- systematic_resampling_n(normweights, nparticles, runif(1))
  ancestors_all <- matrix(nrow = datalength+1, ncol = nparticles)
  ancestors_all[1,] <- ancestors
  # step t >= 1
  for (time in 1:datalength){
    # resampling at every step
    xparticles <- xparticles[,ancestors,drop=FALSE]
    # transition
    xparticles <- rtransition_twisted(xparticles, theta, time, abc)
    xparticles_all[time+1,,] <- xparticles
    # log weights
    logweights <- potential_twisted(xparticles, theta, time, observations[time,], abc)
    logweights_all[time+1,] <- logweights
    maxlogweights <- max(logweights)
    weights <- exp(logweights - maxlogweights)
    # update log likelihood estimator
    ll <- ll + maxlogweights + log(mean(weights))
    ll_all[time+1] <- ll
    # normalize weights
    normweights <- weights / sum(weights)
    # draw ancestors
    ancestors <- systematic_resampling_n(normweights, nparticles, runif(1))
    ancestors_all[time+1,] <- ancestors
  }
  # return all likelihood estimators, particles, weights, ancestors
  return(list(ll = ll, ll_all = ll_all, logweights_all = logweights_all, xparticles_all = xparticles_all, ancestors_all = ancestors_all))
}

## function to update a_t, b_t, c_t at all times, given previous values, and the result of a twisted particle filter
update_abc <- function(pf_results, theta, abc){
  datalength <- dim(pf_results$logweights_all)[1] - 1
  nparticles <- dim(pf_results$logweights_all)[2]
  alpha <- theta[1]
  # sigma2 <- theta[2]
  sigma <- sqrt(sigma2)
  # Now approximate dynamic programming to update a,b,c
  # Initialization, time T+1
  abc_new <- matrix(nrow = datalength+1, ncol = 3)
  time <- datalength
  while (time >= 0){
    x <- pf_results$xparticles_all[time+1,,]
    logMTp1 <- rep(0, nparticles)
    if (time < datalength){
      # M(hat{phi} times psi)
      a_ <- abc[time+2,1] + abc_new[time+2,1]
      b_ <- abc[time+2,2] + abc_new[time+2,2]
      c_ <- abc[time+2,3] + abc_new[time+2,3]
      k_ <- 1/(1/sigma2 + 2*a_)
      logMTp1 <- 0.5*log(k_) - log(sigma) + (0.5 * k_ * (alpha/sigma2 * x - b_)^2 - 0.5 * (alpha/sigma*x)^2 - c_)
      # M(psi)
      a_ <- abc[time+2,1]
      b_ <- abc[time+2,2]
      c_ <- abc[time+2,3]
      k_ <- 1/(1/sigma2 + 2*a_)
      logMTp1 <- logMTp1 - (0.5*log(k_) - log(sigma) + (0.5 * k_ * (alpha/sigma2 * x - b_)^2 - 0.5 * (alpha/sigma*x)^2 - c_))
    }
    logxi <- logMTp1 + pf_results$logweights_all[time+1,]
    x2 <- x^2
    # quadregressioncoef <- as.numeric(lm(-logxi ~ 1 + x + x2)$coef)
    quadregressioncoef <- (.lm.fit(cbind(1,x,x2), -logxi))$coef # coef(.lm.fit(cbind(1,x,x2), -logxi))
    abc_new[time+1,] <- rev(quadregressioncoef) # put the coefficients in the right order
    time <- time - 1
  }
  return(abc + abc_new)
}

## Controlled SMC consists in running twisted particle filters and updating a_t, b_t, c_t
## a number of times (e.g. 3 or 4 times)
csmc <- function(nparticles, theta, observations, niter){
  abc <- matrix(0, nrow = datalength+1, ncol = 3)
  if (niter > 0){
    for (iter in 1:niter){
      pf_results <- particle_filter_twisted(nparticles, theta, observations, abc)
      abc <- update_abc(pf_results, theta, abc)
    }
  }
  pf_results <- particle_filter_twisted(nparticles, theta, observations, abc)
  return(list(ll = pf_results$ll, ll_all = pf_results$ll_all))
}

## prior on theta
## For alpha in (0,1) we put a Uniform prior,
## and for sigma2 we put an inverse gamma prior with parameters (1, 0.1)

## log pdf of inverse Gamma distribution
digamma <- function(x, alpha, beta){
  return(alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x)
}

## log pdf of prior distribution on theta
dprior <- function(theta){
  if (theta[1] < 0 || theta[1] > 1){
    return(-Inf)
  } else {
    return(0)
    # if (theta[2] < 0){
    #   return(-Inf)
    # } else {
    #   return(digamma(theta[2], 1, 0.1))
    # }
  }
}

## with this many particles
nparticles <- 2^8
## and this many iterations per CSMC step
niter <- 3
## and this standard deviation for the random walk proposal
# sd_proposal <- 5 * c(0.002, 0.01) ### values employed for the full data

# stepsizes <- c(0.002, 0.01)
stepsizes <- 0.01

## Markov kernel of the PMMH chain
single_kernel <- function(state){
  position <- state$position
  current_target <- state$current_target
  ## proposal
  sd_proposal <- runif(1, min = 0.1, max = 20) * stepsizes
  # Sigma_chol <- diag(sd_proposal)
  # Sigma_chol_inv <- diag(1/sd_proposal)
  
  proposal <- position + rnorm(1, mean = 0, sd = sd_proposal)
  ## evaluate prior log pdf
  target_proposal <- dprior(proposal)
  if (is.finite(target_proposal)){
    ## only run particle filter if parameter is in the acceptable range
    resll <- try(csmc(nparticles, proposal, observations, niter)$ll)
    if (inherits(resll, "try-error")){
      resll <- -Inf
    }
    if (is.nan(resll)){
      resll <- -Inf
    }
    target_proposal <- target_proposal + resll
  }
  ## accept or not
  if (log(runif(1)) < (target_proposal - current_target)){
    return(list(position = proposal, current_target = target_proposal,
                accept = TRUE))
  } else {
    return(list(position = position, current_target = current_target,
                accept = FALSE))
  }
}

## Markov kernel of the coupled PMMH chain
coupled_kernel <- function(state1, state2){
  position1 <- state1$position
  current_target1 <- state1$current_target
  position2 <- state2$position
  current_target2 <- state2$current_target
  sd_proposal <- runif(1, min = 0.1, max = 20) * stepsizes
  # Sigma_chol <- diag(sd_proposal)
  # Sigma_chol_inv <- diag(1/sd_proposal)
  proposal_value <- unbiasedpoisson::rnorm_reflectionmax(position1, position2, sd_proposal)
  proposal1 <- proposal_value$xy[1]
  proposal2 <- proposal_value$xy[2]
  # proposal_value <- unbiasedpoisson::rmvnorm_reflectionmax(position1, position2, Sigma_chol, Sigma_chol_inv)
  # proposal1 <- proposal_value$xy[,1]
  # proposal2 <- proposal_value$xy[,2]
  logU <- log(runif(1))
  if (proposal_value$identical){
    target_proposal1 <- dprior(proposal1)
    if (is.finite(target_proposal1)){
      resll <- try(csmc(nparticles, proposal1, observations, niter)$ll)
      if (inherits(resll, "try-error")){
        resll <- -Inf
      }
      if (is.nan(resll)){
        resll <- -Inf
      }
      target_proposal1 <- resll + target_proposal1
    }
    target_proposal2 <- target_proposal1
  } else { ## proposals are different
    target_proposal1 <- dprior(proposal1)
    if (is.finite(target_proposal1)){
      resll <- try(csmc(nparticles, proposal1, observations, niter)$ll)
      if (inherits(resll, "try-error")){
        resll <- -Inf
      }
      if (is.nan(resll)){
        resll <- -Inf
      }
      target_proposal1 <- target_proposal1 + resll
    }
    target_proposal2 <- dprior(proposal2)
    if (is.finite(target_proposal2)){
      resll <- try(csmc(nparticles, proposal2, observations, niter)$ll)
      if (inherits(resll, "try-error")){
        resll <- -Inf
      }
      if (is.nan(resll)){
        resll <- -Inf
      }
      target_proposal2 <- target_proposal2 + resll
    }
  }
  accept1 <- FALSE; accept2 <- FALSE
  if (is.finite(target_proposal1)){
    accept1 <- (logU < (target_proposal1 - current_target1))
  }
  if (is.finite(target_proposal2)){
    accept2 <- (logU < (target_proposal2 - current_target2))
  }
  if (accept1){
    position1 <- proposal1
    current_target1 <- target_proposal1
  }
  if (accept2){
    position2 <- proposal2
    current_target2 <- target_proposal2
  }
  identical_ <- proposal_value$identical && accept1 && accept2
  return(list(state1 = list(position = position1, current_target = current_target1),
              state2 = list(position = position2, current_target = current_target2),
              identical = identical_))
}

rinit <- function(theta){
  if (missing(theta)){
    # uniform on [0,1] for parameter
    theta <- c(runif(1, min = 0, max = 1)) 
    # theta <- c(runif(1, min = 0, max = 1), 
    #            runif(1, min = 0, max = 1)) 
  }
  resll <- try(csmc(nparticles, theta, observations, niter)$ll)
  if (inherits(resll, "try-error")){
    resll <- -Inf
  }
  return(list(position = theta, current_target = resll))
}

## test function
h <- function(x) x

## arbitrary state for fishy function estimation
x_0 <- c(0.5)
state_x_0 <- rinit(x_0)


