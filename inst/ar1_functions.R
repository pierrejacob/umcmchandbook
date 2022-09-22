## AR(1) process
## X_t = phi X_{t-1} + W_t, where W_t ~ Normal(0,1)
## invariant distribution is Normal(mean = 0, variance = 1/(1-phi^2))

phi <- 0.99

## Markov kernel
single_kernel <- function(state){
  return(list(position = phi * state$position + rnorm(1, 0, 1)))
}

## Coupled Markov kernel
## using "reflection-maximal coupling" implemented in the function 'rnorm_reflectionmax'
coupled_kernel <- function(state1, state2){
  nextvalues <- unbiasedpoisson::rnorm_reflectionmax(mu1 = phi * state1$position, 
                                                     mu2 = phi * state2$position, 
                                                     sigma = 1)
  return(list(state1 = list(position = nextvalues$xy[1]), 
              state2 = list(position = nextvalues$xy[2]), 
              identical = nextvalues$identical))
}

## initial distribution Normal(0, 4^2)
rinit <- function(position){
  if (missing(position)){
    position <- rnorm(1, mean = 0, sd = 4)
  }
  return(list(position = position))
}

## test function x \mapsto x
h <- function(x){ return(x) }

## arbitrary point used in the definition of the fishy function
x_0 <- 0
state_x_0 <- rinit(x_0)


# # # ## test
# timehorizon <- 1e5
# history <- rep(0, timehorizon)
# state <- rinit()
# for (timeindex in 1:timehorizon){
#   state <- single_kernel(state)
#   history[timeindex] <- state$position
# }
# ##
# matplot(history, type = 'l')
# acf(history)
# # acf(history)$acf[,,1]
# # 
# hist(history, prob = TRUE, nclass = 50, xlab = 'x', main = '')
# curve(dnorm(x, 0, sd = sqrt(1/(1-phi^2))), add = TRUE, lwd = 1.3)
# 
# set.seed(4)
# lag <- 5
# cchain <- sample_coupled_chains(single_kernel, coupled_kernel, rinit, lag = lag, m = 100)
# cchain$meetingtime
# matplot(cchain$samples1, type = 'l', ylim = c(-20, 20))
# matplot(x = (lag+1):(lag+nrow(cchain$samples2)), y = cchain$samples2, type = 'l', add = T,
#         lty = 2, col = 'red')
# 
# ## histogram of meeting times
# nrep <- 1e3
# lag <- 100
# meetingtime_runs <- foreach (irep = 1:nrep) %dorng% {
#   sample_meetingtime(single_kernel, coupled_kernel, rinit, lag = lag)
# }
# meetingtimes <- sapply(meetingtime_runs, function(x) x$meetingtime)
# qplot(x=meetingtimes-lag, geom = 'histogram')
# 
# niterations <- 200
# ubounds <- sapply(1:niterations, function(t) tv_upper_bound(unlist(meetingtimes), lag, t))
# g_tvbounds <- qplot(x = 1:niterations, y = ubounds, geom = "line") +
#   ylab("TV distance") + xlab("time") + ylim(0,1)
# g_tvbounds
