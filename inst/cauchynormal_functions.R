## This script implements couplings
## for the example described in Example 3.1 of 
## Robert, Convergence control methods for MCMC algorithms, Stat Sci 1996

### define kernels
sigma2 <- 100
n <- 3
xobs <- c(-8, 8, 17)

## targetpdf
unnormalizedlogpdf <- function(theta){
  return(-theta^2/(2*sigma2) - sum(log((1+(theta-xobs)^2))))
}

## single kernel
single_kernel <- function(state){
  theta <- state$position
  etas <- rexp(n, rate = 0.5*(1+(theta-xobs)^2))
  var_ <- 1./(sum(etas) + 1./sigma2)
  newtheta <- rnorm(1, mean = var_ * sum(etas * xobs), sd = sqrt(var_))
  return(list(position = newtheta))
}

## coupled kernel that can trigger meetings
coupled_kernel <- function(state1, state2){
  theta1 <- state1$position
  theta2 <- state2$position
  us <- runif(n)
  etas1 <- -log(us) / (0.5*(1+(theta1-xobs)^2))
  etas2 <- -log(us) / (0.5*(1+(theta2-xobs)^2))
  var1_ <- 1./(sum(etas1) + 1./sigma2)
  var2_ <- 1./(sum(etas2) + 1./sigma2)
  mean1_ <- var1_ * sum(etas1 * xobs)
  mean2_ <- var2_ * sum(etas2 * xobs)
  coupledthetas <- unbiasedpoisson::rnorm_max_coupling(mean1_, mean2_, sqrt(var1_), sqrt(var2_))
  return(list(state1 = list(position = coupledthetas$xy[1]), 
              state2 = list(position = coupledthetas$xy[2]),
              identical = coupledthetas$identical))
}

## initial distribution
rinit <- function(x){
  if (missing(x)){
    return(list(position = rnorm(1, 0, 1)))
  } else {
    return(list(position = x))
  }
}


## test function
h <- function(x){
  return(x)
}

## arbitrary point used in the definition of the fishy function
x_0 <- 0
state_x_0 <- rinit(x_0)

## MRTH implementation
target <- function(x){
  if (is.null(dim(x))){
    x <- matrix(x, 1, 1)
  }
  unnormalizedlogpdf(x[1,1])
} 
rinit_mrth <- function(x){
  if (missing(x)){
    x <- rnorm(1, 0, 1)
    return(list(position = matrix(x, 1, 1), current_pdf = target(matrix(x, 1, 1))))
  } else {
    return(list(position = matrix(x, 1, 1), current_pdf = target(matrix(x, 1, 1))))
  }
}
ks <- get_mrth_kernels(target, Sigma_proposal = 10^2)
single_kernel_mrth <- ks$single_kernel
coupled_kernel_mrth <- ks$coupled_kernel

state_x_0_mrth <- rinit_mrth(x_0)

# 
## test
# timehorizon <- 1e5
# history <- rep(0, timehorizon)
# state <- rinit()
# for (timeindex in 1:timehorizon){
#   state <- single_kernel(state)
#   history[timeindex] <- state$position
# }
# #
# matplot(history[1:1e3], type = 'l')
# hist(history, prob = TRUE, nclass = 200, xlab = 'x', main = '', ylim = c(0,0.3))
# # 
# # 
# history_mrth <- rep(0, timehorizon)
# state <- rinit_mrth()
# naccepts <- 0
# for (timeindex in 1:timehorizon){
#   state <- single_kernel_mrth(state)
#   naccepts <- naccepts + state$accept
#   history_mrth[timeindex] <- state$position
# }
# cat("acceptance rate:", naccepts/timehorizon*100, "%\n")
# #
# matplot(history_mrth[1:1e3], type = 'l')
# 
# 
# hist(history, prob = TRUE, nclass = 200, xlab = 'x', main = '', ylim = c(0,0.3), col = rgb(1,1,0,0.5))
# hist(history_mrth, add = T, prob = TRUE, nclass = 200, col = rgb(1,0,0,0.5))
# 
# 
# # targetpdf <- function(z) sapply(z, function(t) exp(unnormalizedlogpdf(t)))
# # normalizingconstant <- integrate(f = targetpdf,lower = -20, upper = 40)$val
# # curve(targetpdf(x)/normalizingconstant, add = TRUE, lwd = 1.3, n = 1000)
# 
# # set.seed(6)
# # lag <- 5
# # cchain <- sample_coupled_chains(single_kernel, coupled_kernel, rinit, lag = lag, m = 100)
# # cchain$meetingtime
# # matplot(cchain$samples1, type = 'l', ylim = c(-20, 20))
# # matplot(x = (lag+1):(lag+nrow(cchain$samples2)), y = cchain$samples2, type = 'l', add = T,
# #         lty = 2, col = 'red')
# # 
# histogram of meeting times
# nrep <- 1e3
# lag <- 100
# meetingtime_runs <- foreach (irep = 1:nrep) %dorng% {
#   sample_meetingtime(single_kernel, coupled_kernel, rinit, lag = lag)
# }
# meetingtimes <- sapply(meetingtime_runs, function(x) x$meetingtime)
# qplot(x=meetingtimes-lag, geom = 'histogram')
# 
# meetingtime_runs_mrth <- foreach (irep = 1:nrep) %dorng% {
#   sample_meetingtime(single_kernel_mrth, coupled_kernel_mrth, rinit_mrth, lag = lag)
# }
# meetingtimes_mrth <- sapply(meetingtime_runs_mrth, function(x) x$meetingtime)
# qplot(x=meetingtimes_mrth-lag, geom = 'histogram')
# 
# # # # 
# niterations <- 100
# ubounds <- sapply(1:niterations, function(t) tv_upper_bound(unlist(meetingtimes), lag, t))
# ubounds_mrth <- sapply(1:niterations, function(t) tv_upper_bound(unlist(meetingtimes_mrth), lag, t))
# g_tvbounds <- qplot(x = 1:niterations, y = ubounds, geom = "line") +
#   ylab("TV distance") + xlab("time") + ylim(0,1)
# g_tvbounds + geom_line(aes(y = ubounds_mrth), linetype = 2)
