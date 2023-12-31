## implements reservoir sampling
#'@export
sample_coupled_chains_and_fish <- function(single_kernel, coupled_kernel, rinit, 
                                                       h, k = 0, m = 1, lag = 1, max_iterations = Inf,
                                                       x_0 = NULL, natoms = 1){
  tictoc::tic("coupled chains and fish")
  if (k > m){
    print("error: k has to be less than m")
    return(NULL)
  }
  if (lag > m){
    print("error: lag has to be less than m")
    return(NULL)
  }
  # starttime <- Sys.time()
  ## format x_0
  state_x_0 <- rinit(x_0)
  dimh <- length(h(state_x_0$position))
  ##### generate first signed measure
  indexselectedatoms <- rep(NA, natoms)
  selectedatoms <- list()
  atomcounter <- 0 # counts atom in the signed measure
  h_at_atoms <- matrix(NA, nrow = dimh, ncol = natoms)
  htilde_at_atoms <- matrix(NA, nrow = dimh, ncol = natoms)
  weight_at_atoms <- rep(NA, natoms)
  ##
  state1 <- rinit(); state2 <- rinit()
  ## initialize selected atoms
  for (iatom in 1:natoms){
    selectedatoms[[iatom]] <- state1
  }
  # mcmcestimator computes the sum of h(X_t) for t=k,...,m
  mcmcestimator <- h(state1$position)
  mcmcestimator_h2 <- h(state1$position)^2
  if (k > 0){
    mcmcestimator <- rep(0, dimh)
    mcmcestimator_h2 <- rep(0, dimh)
    for (iatom in 1:natoms){ selectedatoms[[iatom]] <- NULL }
  } else {
    ## selected atoms might correspond to state1
    atomcounter <- atomcounter + 1
    for (iatom in 1:natoms){
      selectedatoms[[iatom]] <- state1
      indexselectedatoms[iatom] <- atomcounter 
      weight_at_atoms[iatom] <- 1
    }
  }
  # correction computes the sum of w_t * (h(X_{t}) - h(Y_{t-lag})) for t=k+lag,..., tau - 1
  # where w_t = { floor((time-k) / lag) - ceiling(max(lag, time-m)/lag) + 1 } / (m-k+1)
  correction <- rep(0, dimh)
  correction_h2 <- rep(0, dimh)
  time <- 0
  for (t in 1:lag){
    time <- time + 1
    state1 <- single_kernel(state1)
    if (time >= k){
      mcmcestimator <- mcmcestimator + h(state1$position)
      mcmcestimator_h2 <- mcmcestimator_h2 + h(state1$position)^2
      ## selected atoms might correspond to state1
      atomcounter <- atomcounter + 1
      for (iatom in 1:natoms){
        if (runif(1) < 1/atomcounter){ # reservoir sampling mechanism
          indexselectedatoms[iatom] <- atomcounter 
          selectedatoms[[iatom]] <- state1
          weight_at_atoms[iatom] <- 1
        }
      }
    }
  }
  if (time >= k + lag){
    correction <- correction + (floor((time-k) / lag) - ceiling(max(lag, time-m)/lag) + 1) * (h(state1$position) - h(state2$position))
    correction_h2 <- correction_h2 + (floor((time-k) / lag) - ceiling(max(lag, time-m)/lag) + 1) * (h(state1$position)^2 - h(state2$position)^2)
    ## selected atoms might correspond to state1 or state2
    atomcounter <- atomcounter + 1
    for (iatom in 1:natoms){
      if (runif(1) < 1/atomcounter){ # reservoir sampling mechanism
        indexselectedatoms[iatom] <- atomcounter 
        selectedatoms[[iatom]] <- state1
        weight_at_atoms[iatom] <- (floor((time-k) / lag) - ceiling(max(lag, time-m)/lag) + 1)
      }
    }
    atomcounter <- atomcounter + 1
    for (iatom in 1:natoms){
      if (runif(1) < 1/atomcounter){ # reservoir sampling mechanism
        indexselectedatoms[iatom] <- atomcounter 
        selectedatoms[[iatom]] <- state2
        weight_at_atoms[iatom] <- -(floor((time-k) / lag) - ceiling(max(lag, time-m)/lag) + 1)
      }
    }
  }
  meetingtime <- Inf
  # time here is equal to lag; at this point we have X_lag,Y_0 and we are going to generate successively X_{t},Y_{t-lag} where time t is >= lag+1
  while ((time < max(meetingtime, m)) && (time < max_iterations)){
    time <- time + 1 # time is lag+1,lag+2,...
    if (is.finite(meetingtime)){ # chains have already met
      state1 <- single_kernel(state1)
      state2 <- state1
      if (k <= time && time <= m){
        mcmcestimator <- mcmcestimator + h(state1$position)
        mcmcestimator_h2 <- mcmcestimator_h2 + h(state1$position)^2
        atomcounter <- atomcounter + 1
        for (iatom in 1:natoms){
          if (runif(1) < 1/atomcounter){ # reservoir sampling mechanism
            indexselectedatoms[iatom] <- atomcounter 
            selectedatoms[[iatom]] <- state1
            weight_at_atoms[iatom] <- 1
          }
        }
      }
    } else { # chains have not met yet
      res_coupled_kernel <- coupled_kernel(state1, state2)
      state1 <- res_coupled_kernel$state1
      state2 <- res_coupled_kernel$state2
      if (res_coupled_kernel$identical){
        meetingtime <- time
      }
      if (k <= time && time <= m){ # update MCMC estimator 
        mcmcestimator <- mcmcestimator + h(state1$position)
        mcmcestimator_h2 <- mcmcestimator_h2 + h(state1$position)^2
        atomcounter <- atomcounter + 1
        for (iatom in 1:natoms){
          if (runif(1) < 1/atomcounter){ # reservoir sampling mechanism
            indexselectedatoms[iatom] <- atomcounter 
            selectedatoms[[iatom]] <- state1
            weight_at_atoms[iatom] <- 1
          }
        }
        
      }
      if (time >= k + lag){ # update bias correction 
        correction <- correction + (floor((time-k) / lag) - ceiling(max(lag, time-m)/lag) + 1) * (h(state1$position) - h(state2$position))
        correction_h2 <- correction_h2 + (floor((time-k) / lag) - ceiling(max(lag, time-m)/lag) + 1) * (h(state1$position)^2 - h(state2$position)^2)
        ## selected atoms might correspond to state1 or state2
        atomcounter <- atomcounter + 1
        for (iatom in 1:natoms){
          if (runif(1) < 1/atomcounter){
            indexselectedatoms[iatom] <- atomcounter 
            selectedatoms[[iatom]] <- state1
            weight_at_atoms[iatom] <- (floor((time-k) / lag) - ceiling(max(lag, time-m)/lag) + 1)
          }
        }
        atomcounter <- atomcounter + 1
        for (iatom in 1:natoms){
          if (runif(1) < 1/atomcounter){
            indexselectedatoms[iatom] <- atomcounter 
            selectedatoms[[iatom]] <- state2
            weight_at_atoms[iatom] <- -(floor((time-k) / lag) - ceiling(max(lag, time-m)/lag) + 1)
          }
        }
      }
    }
  }
  ## now that we have obtained 'natoms' uniformly selected atoms, 
  ## evaluate h at these atoms and estimate fishy function
  cost_fishyestimation <- rep(0, natoms)
  for (iatom in 1:natoms){
    h_at_atoms[,iatom] <- h(selectedatoms[[iatom]]$position)
    upf <- sample_unbiasedfishy(coupled_kernel, h, selectedatoms[[iatom]], state_x_0, max_iterations = max_iterations)
    htilde_at_atoms[,iatom] <- upf$estimator
    cost_fishyestimation[iatom] <- 2 * upf$meetingtime
  }
  ## finalize unbiased estimator of pi(h) and pi(h^2)
  uestimator <- mcmcestimator + correction
  uestimator_h2 <- mcmcestimator_h2 + correction_h2
  ## cost of obtaining signed measure 
  costsignedmeasure <- lag + 2*(meetingtime - lag) + max(0, time - meetingtime)
  elapsed <- tictoc::toc(quiet = T)
  tictoc::tic.clear()
  return(list(uestimator = uestimator / (m - k + 1),
              uestimator_h2 = uestimator_h2 / (m - k + 1),
              indexselectedatoms = indexselectedatoms,
              selectedatoms = selectedatoms,
              h_at_atoms = h_at_atoms,
              htilde_at_atoms = htilde_at_atoms,
              weight_at_atoms = weight_at_atoms / (m - k + 1),
              atomcounter = atomcounter,
              meetingtime = meetingtime, 
              costsignedmeasure = costsignedmeasure, 
              cost_fishyestimation = cost_fishyestimation, ## cost of _each_ fishy estimator
              elapsedtime = as.numeric(elapsed$toc-elapsed$tic)
        ))
}
