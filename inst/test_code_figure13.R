##########
## Code as in the paper
##########

# Metropolis-Rosenbluth-Teller-Hastings transition with Normal proposals
mrth = function(x, U, sigma)
{
  # proposal = current location + Normal(0,sigma^2)
  xprop = x + sigma * rnorm(length(x))
  # log Uniform to accept/reject proposals
  logu = log(runif(1))
  # return state according to decision to accept or not
  return(if (logu < (U(x) - U(xprop))) xprop else x)
}
# coupling of MRTH transition with Normal proposals
coupledmrth = function(x, y, U, sigma)
{
  # draw proposals using maximal coupling of Bou-Rabee, Eberle and Zimmer (AAP 2020)
  xstd = rnorm(length(x)) # standard Normal variables
  z = (x - y) / sigma # length(sigma) could be 1 or length(x)
  e = z / sqrt(sum(z^2)) # normalise
  logu = log(runif(1)) 
  sameprop = (logu < sum(dnorm(xstd + z, log = TRUE) - dnorm(xstd, log = TRUE)))
  ystd = if (sameprop) xstd + z else xstd - 2 * sum(e * xstd) * e
  # xprop is marginally Normal(x,sigma^2) 
  xprop = x + sigma * xstd
  # yprop is marginally Normal(y,sigma^2)
  yprop = y + sigma * ystd
  # log Uniform to accept/reject proposals
  logu = log(runif(1))
  # decision to accept or not
  xaccept = (logu < (U(x) - U(xprop)))
  yaccept = (logu < (U(y) - U(yprop)))
  # return state according to decision
  return(list(nextx = if (xaccept) xprop else x,
              nexty = if (yaccept) yprop else y,
              nextxequalsnexty = sameprop && xaccept && yaccept))
}



##########
## test 
##########

sigma2 <- 100
## target probability density function (pdf)
xobs <- c(-8, 8, 17)
unnormalizedlogpdf <- function(theta){
  return(-theta^2/(2*sigma2) - sum(log((1+(theta-xobs)^2))))
}
U <- function(x) -unnormalizedlogpdf(x)
##
sigma_prop <- 2
rinit <- function() rnorm(1, -10, 10)

## run MCMC
set.seed(1)
nmcmc <- 50000
xchain <- rep(0, nmcmc)
xchain[1] <- rinit()
for (imcmc in 2:nmcmc){
  xchain[imcmc] <- mrth(xchain[imcmc-1], U, sigma_prop)
}

## trace plot
plot(xchain[1:500], type = 'l')

## marginal density post-burnin vs target density
targetpdf <- function(z) sapply(z, function(t) exp(unnormalizedlogpdf(t)))
normalizingconstant <- integrate(f = targetpdf,lower = -30, upper = 40)$val
xseq <- seq(from = -30, to = 40, length.out = 1000)
plot(xseq, sapply(xseq, function(z) targetpdf(z)/normalizingconstant), type = 'l', xlab = 'x', ylab = 'density', col = 'red')
lines(density(xchain[500:nmcmc]), lty = 2)

## generate meeting times
lag <- 100
nrep <- 1000
meetingtimes <- rep(0, nrep)
for (irep in 1:nrep){
  x <- rinit()
  y <- rinit()
  for (t in 1:lag) x <- mrth(x, U, sigma_prop)
  time <- lag
  while (time < 1e5){
    time <- time + 1
    res <- coupledmrth(x, y, U, sigma_prop)
    x <- res$nextx
    y <- res$nexty
    if (res$nextxequalsnexty) break 
  }
  meetingtimes[irep] <- time
}
## histogram of meeting times - lag
hist(meetingtimes-lag, xlab = "meeting times - lag", nclass = 50)

## TV upper bounds
plot(sapply(1:max(meetingtimes-lag), function(t) return(mean(pmax(0, ceiling((meetingtimes-lag-t) / lag))))), type = 'l', xlab = "Iteration", ylab = "TV upper bounds (log scale)", log = 'y')


##########
## two dimensional example
##########

U <- function(x) (x[1]-1)^2 / (2*1) + (x[2]+1)^2 / (2*2)
##
sigma_prop <- c(2.5, 1.4)
rinit <- function() rnorm(2, -10, 10)

## run MCMC
set.seed(1)
nmcmc <- 50000
xchain <- matrix(0, nrow = nmcmc, ncol = 2)
xchain[1,] <- rinit()
for (imcmc in 2:nmcmc){
  xchain[imcmc,] <- mrth(xchain[imcmc-1,], U, sigma_prop)
}

## trace plot
matplot(xchain[1:500,], type = 'l')

## marginal density post-burnin vs target density
hist(xchain[500:nmcmc,1], prob = TRUE, nclass = 100)
curve(dnorm(x, mean = 1, sd = 1), add = TRUE)
hist(xchain[500:nmcmc,2], prob = TRUE, nclass = 100)
curve(dnorm(x, mean =  -1, sd = sqrt(2)), add = TRUE)

## generate meeting times
lag <- 50
nrep <- 1000
meetingtimes <- rep(0, nrep)
for (irep in 1:nrep){
  x <- rinit()
  y <- rinit()
  for (t in 1:lag) x <- mrth(x, U, sigma_prop)
  time <- lag
  while (time < 1e5){
    time <- time + 1
    res <- coupledmrth(x, y, U, sigma_prop)
    x <- res$nextx
    y <- res$nexty
    if (res$nextxequalsnexty) break 
  }
  meetingtimes[irep] <- time
}
## histogram of meeting times - lag
hist(meetingtimes-lag, xlab = "meeting times - lag", nclass = 50)
## TV upper bounds
plot(sapply(1:(max(meetingtimes-lag)), function(t) return(mean(pmax(0, ceiling((meetingtimes-lag-t) / lag))))), type = 'l', xlab = "Iteration", ylab = "TV upper bounds (log scale)", log = 'y')


