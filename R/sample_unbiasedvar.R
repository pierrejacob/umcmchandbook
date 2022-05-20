## sample unbiased estimator of asymptotic variance
#'@rdname sample_unbiasedvar
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
#' \item pih:
#' \item cost_fishyestimation:
#' \item cost_umcmc:
#' \item cost:
#'}
#'@export
sample_unbiasedvar <- function(single_kernel, coupled_kernel, rinit, h, 
                          k = 0, m = 1, lag = 1, x_0 = NULL, natoms = 1){
  tictoc::tic("unbiased asymptotic variance estimator")
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
  xis <- rep(1/length(cc1_df$weight),length(cc1_df$weight))
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
  xis <- rep(1/length(cc2_df$weight),length(cc2_df$weight))
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
  ## estimator of the "covariance" term XXXX
  ## estimator of the variance of h under pi
  ## estimator of pi(h)
  ## and cost measured in number of iterations
  elapsed <- tictoc::toc(quiet = T)
  tictoc::tic.clear()
  return(list(estimator = asymvar_estimator, covariance_term = covariance_term,
              varh = varh, pih = pi_h_average, cost_fishyestimation = cost_fishyestimation,
              cost_umcmc = cost_umcmc, cost = cost_umcmc + cost_fishyestimation,
              elapsedtime = as.numeric(elapsed$toc-elapsed$tic)))
}