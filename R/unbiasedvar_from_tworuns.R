
## assemble two runs of "sample_coupled_chains_and_fish" into
## an unbiased estimator of v(P,h)
#'@export
unbiasedvar_from_tworuns <- function(run1, run2, natoms = 1){
  tictoc::tic("unbiased asymptotic variance estimator")
  ## if natoms is a vector, run algorithm with R=max(natoms)
  natoms <- sort(natoms)
  natoms_max <- natoms[length(natoms)]
  ## estimator of pi(h^2) as average of two estimators
  pi_hsquared_ <- (1/2)*(run1$uestimator_h2 + run2$uestimator_h2)
  ## estimator of pi(h)^2 as product of two independent estimators of pi(h)
  pih_deepbreath_squared_ <- run1$uestimator * run2$uestimator
  ## estimate of variance under pi
  varh <- pi_hsquared_ - pih_deepbreath_squared_
  ## term involving fishy function estimates
  ## computed for each 'natoms' if 'natoms' is a vector
  fishyterm1 <- run1$atomcounter * t(apply(X = (run1$h_at_atoms - run2$uestimator) * run1$htilde_at_atoms, MARGIN = 1, FUN = function(v) cumsum(run1$weight_at_atoms * v)))
  fishyterm2 <- run2$atomcounter * t(apply(X = (run2$h_at_atoms - run1$uestimator) * run2$htilde_at_atoms, MARGIN = 1, FUN = function(v) cumsum(run2$weight_at_atoms * v)))
  fishyterms <- t(t(fishyterm1 + fishyterm2)/(1:natoms_max))[,natoms,drop=F]
  estimator <- - varh + fishyterms ## compute unbiased asymptotic variance estimator
  cost_fishyterms <- cumsum(run1$cost_fishyestimation + run2$cost_fishyestimation)[natoms]
  cost <- run1$costsignedmeasure + run2$costsignedmeasure + cost_fishyterms
  elapsed <- tictoc::toc(quiet = T)
  tictoc::tic.clear()
  return(list(estimator = estimator, fishyterms = fishyterms,
              varh = varh, pih = cbind(run1$uestimator, run2$uestimator), 
              cost_fishyterms = cost_fishyterms, 
              cost_umcmc = run1$costsignedmeasure + run2$costsignedmeasure,
              cost = cost, 
              elapsedtime = as.numeric(elapsed$toc-elapsed$tic)))
}
