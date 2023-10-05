## sample unbiased estimator of asymptotic variance

### WARNING this works only if a chain can "rinit" again given the position
### probably better to focus on sample_unbiasedvar_reservoir which does not have this issue

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
#' \item fishyterms:
#' \item varh:
#' \item pih:
#' \item cost_fishyterms:
#' \item cost_umcmc:
#' \item cost:
#'}
#'@export
sample_unbiasedvar <- function(single_kernel, coupled_kernel, rinit, h, 
                          k = 0, m = 1, lag = 1, x_0 = NULL, natoms = 1, max_iterations = Inf){
  tictoc::tic("unbiased asymptotic variance estimator")
  ## run two coupled chains
  cc1_ <-   umcmchandbook::sample_coupled_chains(single_kernel, coupled_kernel, rinit, m = m, lag = lag, max_iterations = max_iterations)
  cc2_ <-   umcmchandbook::sample_coupled_chains(single_kernel, coupled_kernel, rinit, m = m, lag = lag, max_iterations = max_iterations)
  cost_umcmc <- cc1_$cost + cc2_$cost
  ## from both signed measures, get unbiased estimator of pi(h) and pi(h^2)
  pi_h1_        <- umcmchandbook::H_bar(cc1_, h = h, k = k, m = m)
  pi_h2_        <- umcmchandbook::H_bar(cc2_, h = h, k = k, m = m)
  pi_hsquared1_ <- umcmchandbook::H_bar(cc1_, h = function(x) h(x)^2, k = k, m = m)
  pi_hsquared2_ <- umcmchandbook::H_bar(cc2_, h = function(x) h(x)^2, k = k, m = m)
  dimh <- length(pi_h1_)
  ## estimator of pi(h^2) as average of two estimators
  pi_hsquared_ <- (1/2)*(pi_hsquared1_ + pi_hsquared2_)
  ## estimator of pi(h)^2 as product of two independent estimators of pi(h)
  pih_deepbreath_squared_ <- pi_h1_ * pi_h2_
  ## estimate of variance under pi
  varh <- pi_hsquared_ - pih_deepbreath_squared_
  fishyterms <- rep(0, dimh) 
  cost_fishyterms <- 0
  ## next get fishy function estimator for atoms in cc1_
  ## if natoms is a vector, run algorithm with R=max(natoms)
  natoms <- sort(natoms)
  natoms_max <- natoms[length(natoms)]
  ## convert to data frame and prune identical atoms 
  cc1_df <- umcmchandbook:::c_chains_to_dataframe(cc1_, k, m, dopar = FALSE, prune = TRUE)
  natomsincc1 <- length(cc1_df$weight)
  xis <- rep(1./natomsincc1, natomsincc1)
  h_at_atoms1 <- matrix(NA, nrow = dimh, ncol = natoms_max)
  htilde_at_atoms1 <- matrix(NA, nrow = dimh, ncol = natoms_max)
  weight_at_atoms1 <- rep(NA, natoms_max)
  cost_fishyestimation1 <- rep(0, natoms_max)
  for (iatom in 1:natoms_max){
    ## select an atom at random
    index_ell <- sample(x = 1:natomsincc1, size = 1, prob = xis)
    atom_ell <- as.numeric(cc1_df[index_ell,4:dim(cc1_df)[2]])
    ## start two chains
    state_atom_ <- rinit(as.numeric(atom_ell))
    state_x_0 <- rinit(x_0)
    ## estimate fishy function at that atom
    upfrun <- sample_unbiasedfishy(coupled_kernel, h, state_atom_, state_x_0, max_iterations = max_iterations)
    h_at_atoms1[,iatom] <- h(atom_ell)
    htilde_at_atoms1[,iatom] <- upfrun$estimator
    cost_fishyestimation1[iatom] <- 2 * upfrun$meetingtime
    weight_at_atoms1[iatom] <- cc1_df$weight[index_ell]
  }
  ## next get fishy function estimator for atoms in cc2_
  ## convert to data frame and prune identical atoms 
  cc2_df <- umcmchandbook:::c_chains_to_dataframe(cc2_, k, m, dopar = FALSE, prune = TRUE)
  natomsincc2 <- length(cc2_df$weight)
  xis <- rep(1./natomsincc2, natomsincc2)
  h_at_atoms2 <- matrix(NA, nrow = dimh, ncol = natoms_max)
  htilde_at_atoms2 <- matrix(NA, nrow = dimh, ncol = natoms_max)
  weight_at_atoms2 <- rep(NA, natoms_max)
  cost_fishyestimation2 <- rep(0, natoms_max)
  for (iatom in 1:natoms_max){
    ## select an atom at random
    index_ell <- sample(x = 1:natomsincc2, size = 1, prob = xis)
    atom_ell <- as.numeric(cc2_df[index_ell,4:dim(cc2_df)[2]])
    ## start two chains
    state_atom_ <- rinit(as.numeric(atom_ell))
    state_x_0 <- rinit(x_0)
    ## estimate fishy function at that atom
    upfrun <- sample_unbiasedfishy(coupled_kernel, h, state_atom_, state_x_0, max_iterations = max_iterations)
    h_at_atoms2[,iatom] <- h(atom_ell)
    htilde_at_atoms2[,iatom] <- upfrun$estimator
    cost_fishyestimation2[iatom] <- 2 * upfrun$meetingtime
    weight_at_atoms2[iatom] <- cc2_df$weight[index_ell]
  }
  fishyterm1 <- natomsincc1 * t(apply(X = (h_at_atoms1 - pi_h2_) * htilde_at_atoms1, MARGIN = 1, FUN = function(v) cumsum(weight_at_atoms1 * v)))
  fishyterm2 <- natomsincc2 * t(apply(X = (h_at_atoms2 - pi_h1_) * htilde_at_atoms2, MARGIN = 1, FUN = function(v) cumsum(weight_at_atoms2 * v)))
  fishyterms <- t(t(fishyterm1 + fishyterm2)/(1:natoms_max))[,natoms,drop=F]
  estimator <- - varh + fishyterms ## compute unbiased asymptotic variance estimator
  cost_fishyterms <- cumsum(cost_fishyestimation1 + cost_fishyestimation2)[natoms]
  ## returns estimator of the asymptotic variance
  ## estimator of the "covariance" term with fishy function estimators
  ## estimator of the variance of h under pi
  ## estimator of pi(h)
  ## and cost measured in number of iterations
  elapsed <- tictoc::toc(quiet = T)
  tictoc::tic.clear()
  return(list(estimator = estimator, fishyterms = fishyterms,
              varh = varh, pih = cbind(pi_h1_, pi_h2_), cost_fishyterms = cost_fishyterms,
              cost_umcmc = cost_umcmc, cost = cost_umcmc + cost_fishyterms,
              elapsedtime = as.numeric(elapsed$toc-elapsed$tic)))
}


