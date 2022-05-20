## sample unbiased estimator of asymptotic variance with light memory usage using reservoir sampling
#'@rdname sample_unbiasedvar_reservoir
#'@title Sample unbiased estimator of asymptotic variance in CLT for Markov chain averages
#'@description Sample unbiased estimator of asymptotic variance
#' in the CLT for Markov chain averages.
#'@param single_kernel A list taking a state and returning a state, performing one step of a Markov kernel
#'@param coupled_kernel A list taking two states and returning two states, performing one step of a coupled Markov kernel;
#'it also returns a boolean "identical" indicating whether the two states are identical.
#'@param rinit A list representing the initial state of the chain, that can be given to 'single_kernel'
#'@param h test function h
#'@param k An integer at which to start computing the unbiased estimator
#'@param m A time horizon: the chains are sampled until the maximum between m and the meeting time
#'@param lag A time lag, equal to one by default
#'@param x_0 state fixed arbitrarily to define fishy function (y in the paper)
#'@param natoms number of fishy function evaluations to estimate per signed measure (R in the paper)
#'@return A list with
#'\itemize{
#' \item estimator: 
#' \item covariance_term:
#' \item varh:
#' \item pih: vector with two unbiased estimators of pi(h)
#' \item cost_fishyestimation: cost coming from estimators of fishy function evaluations
#' \item cost: total cost
#' \item elapsedtime: in seconds
#'}
#'@export
sample_unbiasedvar_reservoir <- function(single_kernel, coupled_kernel, rinit, h, 
                                         k = 0, m = 1, lag = 1, x_0 = NULL, natoms = 1){
  tictoc::tic("unbiased asymptotic variance estimator")
  # starttime <- Sys.time()
  ## if natoms is a vector, run algorithm with R=max(natoms)
  natoms <- sort(natoms)
  natoms_max <- natoms[length(natoms)]
  run1 <- unbiasedpoisson:::sample_coupled_chains_and_fish(single_kernel, coupled_kernel, rinit, h = h, k = k, m = m, lag = lag, x_0 = x_0, natoms = natoms_max)
  run2 <- unbiasedpoisson:::sample_coupled_chains_and_fish(single_kernel, coupled_kernel, rinit, h = h, k = k, m = m, lag = lag, x_0 = x_0, natoms = natoms_max)
  ## estimator of pi(h^2) as average of two estimators
  pi_hsquared_ <- (1/2)*(run1$uestimator_h2 + run2$uestimator_h2)
  ## estimator of pi(h)^2 as product of two independent estimators of pi(h)
  pih_deepbreath_squared_ <- run1$uestimator * run2$uestimator
  ## estimate of variance under pi
  varh <- pi_hsquared_ - pih_deepbreath_squared_
  ## term involving fishy function estimates
  ## computed for each 'natoms' if 'natoms' is a vector
  fishyterm1 <- run1$atomcounter * cumsum(run1$weight_at_atoms * (run1$h_at_atoms - run2$uestimator) * run1$htilde_at_atoms)
  fishyterm2 <- run2$atomcounter * cumsum(run2$weight_at_atoms * (run2$h_at_atoms - run1$uestimator) * run2$htilde_at_atoms)
  fishyterms <- ((fishyterm1 + fishyterm2)/(1:natoms_max))[natoms]
  estimator <- - varh + fishyterms ## compute unbiased asymptotic variance estimator
  cost_fishyterms <- cumsum(run1$cost_fishyestimation + run2$cost_fishyestimation)[natoms]
  cost <- run1$costsignedmeasure + run2$costsignedmeasure + cost_fishyterms
  # currentime <- Sys.time()
  # elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currentime) - lubridate::ymd_hms(starttime)), "seconds")
  elapsed <- tictoc::toc(quiet = T)
  tictoc::tic.clear()
  return(list(estimator = estimator, fishyterms = fishyterms,
              varh = varh, pih = c(run1$uestimator, run2$uestimator), 
              cost_fishyterms = cost_fishyterms, cost = cost, 
              elapsedtime = as.numeric(elapsed$toc-elapsed$tic)))
}

