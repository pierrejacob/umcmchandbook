### TAKEN FROM https://github.com/niloyb/CoupledHalfT
### and slightly modified

## Couplings for xi update given eta ##
# Coupled Metropolis-Hastings update of xi given eta
#' xi_update
#' @description Coupled Metropolis-Hastings update of xi given eta
#' @param current_xi_1,current_xi_2 current xi values (positive scalars)
#' @param eta_1,eta_2 current eta values (vector length p)
#' @param X n by p matrix
#' @param X_transpose Pre-calculated transpose of X
#' @param y length n vector
#' @param a0 positive scalar
#' @param b0 positive scalar
#' @param std_MH standard deviation of log-normal MH proposal
#' @return Coupled Metropolis-Hastings update of xi given eta
maxcoupling_xi_update <- function(current_xi_1, eta_1, current_xi_2, eta_2,
                                X, X_transpose, y, a0, b0, std_MH){
  standard_normal <- rnorm(1, mean = 0, sd = 1)
  log_proposed_xi_1 <- standard_normal*sqrt(std_MH) + log(current_xi_1)
  # using max coupling to get the proposal in the MRTH algorithm
  if (dnorm(log_proposed_xi_1, mean = log(current_xi_1), sd = std_MH, log = TRUE) +
      log(runif(1)) < dnorm(log_proposed_xi_1, mean = log(current_xi_2), sd = std_MH, log = TRUE)){
    log_proposed_xi_2 <- log_proposed_xi_1
  } else {
    reject <- TRUE
    y_proposal <- NA
    attempts <- 0
    while(reject){
      attempts <- attempts + 1
      y_proposal <- rnorm(1, mean = log(current_xi_2), sd = std_MH)
      reject <- (dnorm(y_proposal, mean = log(current_xi_2), sd = std_MH, log = TRUE) +
                   log(runif(1)) < dnorm(y_proposal, mean = log(current_xi_1), sd = std_MH, log = TRUE))
    }
    log_proposed_xi_2 <- y_proposal
  }
  proposed_xi_1 <- exp(log_proposed_xi_1)
  proposed_xi_2 <- exp(log_proposed_xi_2)
  
  X_eta_tX_matrix_1 <- X_eta_tX(eta_1, X, X_transpose)
  log_ratio_current_ssr_matrixinv_1 <- log_ratio(current_xi_1, eta_1, X_eta_tX_matrix_1, y, a0, b0)
  log_ratio_proposed_ssr_matrixinv_1 <- log_ratio(proposed_xi_1, eta_1, X_eta_tX_matrix_1, y, a0, b0)
  
  X_eta_tX_matrix_2 <- X_eta_tX(eta_2, X, X_transpose)
  log_ratio_current_ssr_matrixinv_2 <- log_ratio(current_xi_2, eta_2, X_eta_tX_matrix_2, y, a0, b0)
  log_ratio_proposed_ssr_matrixinv_2 <- log_ratio(proposed_xi_2, eta_2, X_eta_tX_matrix_2, y, a0, b0)
  
  log_u <- log(runif(1))
  
  log_accept_prob_1 <- (log_ratio_proposed_ssr_matrixinv_1$log_likelihood - log_ratio_current_ssr_matrixinv_1$log_likelihood) + (log(proposed_xi_1) - log(current_xi_1))
  if (log_u < log_accept_prob_1){
    current_xi_1 <- proposed_xi_1
    log_ratio_current_ssr_matrixinv_1 <- log_ratio_proposed_ssr_matrixinv_1
  }
  
  log_accept_prob_2 <- (log_ratio_proposed_ssr_matrixinv_2$log_likelihood - log_ratio_current_ssr_matrixinv_2$log_likelihood) + (log(proposed_xi_2) - log(current_xi_2))
  if (log_u < log_accept_prob_2){
    current_xi_2 <- proposed_xi_2
    log_ratio_current_ssr_matrixinv_2 <- log_ratio_proposed_ssr_matrixinv_2
  }
  return(list('xi_values' = c(current_xi_1, current_xi_2),
              'log_ratio_ssr_matrix_inv_1' = log_ratio_current_ssr_matrixinv_1,
              'log_ratio_ssr_matrix_inv_2' = log_ratio_current_ssr_matrixinv_2))
}

## Couplings for sigma^2 update given eta ##
digamma <- function(x, alpha, beta){
  return(alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x)
}

rigamma <- function(n, alpha, beta){
  return(1/rgamma(n = n, shape = alpha, rate = beta))
}

maxcoupling_rigamma <- function(alpha1, alpha2, beta1, beta2){
  x <- rigamma(1, alpha1, beta1)
  if (digamma(x, alpha1, beta1) + log(runif(1)) < digamma(x, alpha2, beta2)){
    return(c(x,x))
  } else {
    reject <- TRUE
    y <- NA
    while (reject){
      y <- rigamma(1, alpha2, beta2)
      reject <- (digamma(y, alpha2, beta2) + log(runif(1)) < digamma(y, alpha1, beta1))
    }
    return(c(x,y))
  }
}

# Coupled update of sigma2 given eta
#' xi_update
#' @description Coupled update of sigma2 given eta
#' @param xi_1,xi_2 current xi values (positive scalars)
#' @param eta_1,eta_2 current eta values (vector length p)
#' @param n number of data points
#' @param ssr_1,ssr_2 postiive scalars
#' @param a0 positive scalar
#' @param b0 positive scalar
#' @return Coupled update of sigma2 given eta
maxcoupling_sigma2_update <- function(xi_1, eta_1, xi_2, eta_2, n, ssr_1, ssr_2, a0, b0){
  output <- maxcoupling_rigamma(((n+a0)/2), ((n+a0)/2), (ssr_1/2), (ssr_2/2))
  return(output)
}

## Common random numbers coupling of beta given xi, eta, sigma2 ##
#'@export
crn_beta_update <-
  function(xi_1, sigma2_1, eta_1, xi_2, sigma2_2, eta_2,
           X, X_transpose, y, M_matrix_inverse_1, M_matrix_inverse_2){
  n <- nrow(X)
  p <- nrow(X_transpose)
  # Using same common random numbers for draws on two chains
  random_u <- rnorm(p, 0, 1)
  random_delta <- c(rnorm(n,0,1))
  u_1 = (sqrt(xi_1*eta_1)^-1)*random_u
  v_1 = X%*%u_1 + random_delta
  v_star_1 <- M_matrix_inverse_1%*%(y/sqrt(sigma2_1) - v_1)
  U_1 = (xi_1^-1)*((eta_1^(-1))*(X_transpose))
  u_1 <- u_1 + U_1%*%v_star_1
  beta_parameter_1 <- sqrt(sigma2_1)*(u_1)
  u_2 = (sqrt(xi_2*eta_2)^-1)*random_u
  v_2 = X%*%u_2 + random_delta
  v_star_2 <- M_matrix_inverse_2%*%(y/sqrt(sigma2_2) - v_2)
  U_2 = (xi_2^-1)*((eta_2^(-1))*(X_transpose))
  u_2 <- u_2 + U_2%*%v_star_2
  beta_parameter_2 <- sqrt(sigma2_2)*(u_2)
  return(cbind(beta_parameter_1, beta_parameter_2))
}

#' coupled_half_t_kernel
#' @description Coupled blocked Gibbs kernel for half-t priors
#' @param X n by p matrix
#' @param X_transpose Pre-calculated transpose of X
#' @param y length n vector
#' @param a0 positive scalar
#' @param b0 positive scalar
#' @param std_MH standard deviation of log-normal MH proposal
#' @param xi_1_current,xi_2_current current xi values (positive scalar)
#' @param sigma2_1_current,sigma2_2_current current sigma2 values (positive scalar)
#' @param beta_1_current,beta_2_current current beta values (vector of length p)
#' @param eta_1_current,eta_2_current current eta values (vector of length p)
#' @param nrepeats_eta number of slice sampling steps
#' @param t_dist_df degree of freedom v>=1 for Half-t(v).
#' @return (beta1, beta2, eta1, eta2, sigma2_1, sigma2_2, xi1, xi2, metric_d)
#' @export
coupled_half_t_kernel <- function(X, X_transpose, y, a0=1, b0=1, std_MH=0.8,
                                  xi_1_current, xi_2_current, sigma2_1_current, sigma2_2_current,
                                  beta_1_current, beta_2_current, eta_1_current, eta_2_current,
                                  nrepeats_eta=1,
                                  t_dist_df){
  n <- dim(X)[1]
  p <- dim(X)[2]
  if (nrepeats_eta > 1){
    for (irepeta in 1:(nrepeats_eta-1)){
      eta_sample <- crn_eta_update(xi_1_current, beta_1_current, eta_1_current, sigma2_1_current,
                                   xi_2_current, beta_2_current, eta_2_current, sigma2_2_current, t_dist_df)
      eta_1_current <- eta_sample[,1]
      eta_2_current <- eta_sample[,2]
      
    }
  }
  eta_sample <- swithtocrn_eta_update(xi_1_current, beta_1_current, eta_1_current, sigma2_1_current,
                                      xi_2_current, beta_2_current, eta_2_current, sigma2_2_current, t_dist_df)
  
  eta_1_new <- eta_sample[,1]
  eta_2_new <- eta_sample[,2]
  xi_sample <-
    maxcoupling_xi_update(xi_1_current, eta_1_new, xi_2_current, eta_2_new, 
                          X, X_transpose, y, a0, b0, std_MH)
  xi_1_new <- xi_sample$xi_values[1]
  xi_2_new <- xi_sample$xi_values[2]
  sigma2_sample <-
    maxcoupling_sigma2_update(xi_1_new,eta_1_new,xi_2_new,eta_2_new, n,
                              (xi_sample$log_ratio_ssr_matrix_inv_1)$ssr,
                              (xi_sample$log_ratio_ssr_matrix_inv_2)$ssr, a0, b0)
  sigma2_1_new <- sigma2_sample[1]
  sigma2_2_new <- sigma2_sample[2]
  M_inverse_1 <- xi_sample$log_ratio_ssr_matrix_inv_1$M_matrix_inverse
  M_inverse_2 <- xi_sample$log_ratio_ssr_matrix_inv_2$M_matrix_inverse
  beta <- crn_beta_update(xi_1_new, sigma2_1_new, eta_1_new,
                                  xi_2_new, sigma2_2_new, eta_2_new,
                                  X, X_transpose, y, M_inverse_1, M_inverse_2)
  beta_1_new <- beta[,1]
  beta_2_new <- beta[,2]
  output <- list('beta_1'=beta_1_new, 'beta_2'=beta_2_new,
                 'eta_1'=eta_1_new, 'eta_2'=eta_2_new,
                 'sigma2_1'=sigma2_1_new, 'sigma2_2'=sigma2_2_new,
                 'xi_1'=xi_1_new, 'xi_2'=xi_2_new)
  return(output)
}



#' 
#' ## Full coupled blocked Gibbs samplers ##
#' #' coupled_half_t_kernel
#' #' @description Coupled blocked Gibbs kernel for half-t priors
#' #' @param X n by p matrix
#' #' @param X_transpose Pre-calculated transpose of X
#' #' @param y length n vector
#' #' @param a0 positive scalar
#' #' @param b0 positive scalar
#' #' @param std_MH standard deviation of log-normal MH proposal
#' #' @param xi_1_current,xi_2_current current xi values (positive scalar)
#' #' @param sigma2_1_current,sigma2_2_current current sigma2 values (positive scalar)
#' #' @param beta_1_current,beta_2_current current beta values (vector of length p)
#' #' @param eta_1_current,eta_2_current current eta values (vector of length p)
#' #' @param epsilon_eta eta common random numbers/ maximal coupling coupling threshold
#' #' @param nrepeats_eta number of slice sampling steps
#' #' @param t_dist_df degree of freedom v>=1 for Half-t(v).
#' #' @return (beta1, beta2, eta1, eta2, sigma2_1, sigma2_2, xi1, xi2, metric_d)
#' #' @export
#' coupled_half_t_kernel <- function(X, X_transpose, y, a0=1, b0=1, std_MH=0.8,
#'            xi_1_current, xi_2_current, sigma2_1_current, sigma2_2_current,
#'            beta_1_current, beta_2_current, eta_1_current, eta_2_current,
#'            epsilon_eta = 0.5, nrepeats_eta=1,
#'            t_dist_df, two_scale=TRUE){
#'     n <- dim(X)[1]
#'     p <- dim(X)[2]
#' 
#'     if (two_scale){
#'       # Calculating the metric
#'       # When epsilon_eta >=1, relative_error_delta <= epsilon_eta always.
#'       if (epsilon_eta >= 1){
#'         relative_error_delta <- 1
#'       } else {
#'         relative_error_delta <-
#'           half_t_max_couple_prob(xi_1_current, beta_1_current, eta_1_current, sigma2_1_current,
#'                                  xi_2_current, beta_2_current, eta_2_current, sigma2_2_current,
#'                                  t_dist_df, iterations=1)
#'       }
#'       if (relative_error_delta <= epsilon_eta){ # Using max coupling of 1-step slice sampling when close
#'         eta_sample <-
#'           eta_update_half_t_max_couple(xi_1_current, beta_1_current, eta_1_current, sigma2_1_current,
#'                                        xi_2_current, beta_2_current, eta_2_current, sigma2_2_current,
#'                                        t_dist_df)
#'         # eta_sample <-
#'         #   eta_update_half_t_max_couple_till_you_miss(xi_1_current, beta_1_current, eta_1_current, sigma2_1_current,
#'         #                                              xi_2_current, beta_2_current, eta_2_current, sigma2_2_current,
#'         #                                              t_dist_df)
#'       } else {
#'         for (i in 1:nrepeats_eta){
#'           eta_sample <-
#'             eta_update_half_t_crn_couple(xi_1_current, beta_1_current, eta_1_current, sigma2_1_current,
#'                                          xi_2_current, beta_2_current, eta_2_current, sigma2_2_current,
#'                                          t_dist_df)
#'           eta_1_current <- eta_sample[,1]
#'           eta_2_current <- eta_sample[,2]
#'         }
#'       }
#'     } else {
#'       relative_error_delta <- 1 # Now no need to calculate the metric
#'       eta_sample <-
#'         swithtocrn_eta_update(xi_1_current, beta_1_current, eta_1_current, sigma2_1_current,
#'                                                    xi_2_current, beta_2_current, eta_2_current, sigma2_2_current,
#'                                                    t_dist_df)
#'     }
#'     
#'     eta_1_new <- eta_sample[,1]
#'     eta_2_new <- eta_sample[,2]
#'     xi_sample <-
#'       maxcoupling_xi_update(xi_1_current, eta_1_new, xi_2_current, eta_2_new, 
#'                       X, X_transpose, y, a0, b0, std_MH)
#'     xi_1_new <- xi_sample$xi_values[1]
#'     xi_2_new <- xi_sample$xi_values[2]
#'     sigma2_sample <-
#'       maxcoupling_sigma2_update(xi_1_new,eta_1_new,xi_2_new,eta_2_new,n,
#'                               (xi_sample$log_ratio_ssr_matrix_inv_1)$ssr,
#'                               (xi_sample$log_ratio_ssr_matrix_inv_2)$ssr, a0, b0)
#'     sigma2_1_new <- sigma2_sample[1]
#'     sigma2_2_new <- sigma2_sample[2]
#'     M_inverse_1 <- xi_sample$log_ratio_ssr_matrix_inv_1$M_matrix_inverse
#'     M_inverse_2 <- xi_sample$log_ratio_ssr_matrix_inv_2$M_matrix_inverse
#'     beta <- crn_beta_update(xi_1_new, sigma2_1_new, eta_1_new,
#'                                           xi_2_new, sigma2_2_new, eta_2_new,
#'                                           X, X_transpose, y, M_inverse_1, M_inverse_2)
#'     beta_1_new <- beta[,1]
#'     beta_2_new <- beta[,2]
#'     output <- list('beta_1'=beta_1_new, 'beta_2'=beta_2_new,
#'                    'eta_1'=eta_1_new, 'eta_2'=eta_2_new,
#'                    'sigma2_1'=sigma2_1_new, 'sigma2_2'=sigma2_2_new,
#'                    'xi_1'=xi_1_new, 'xi_2'=xi_2_new,
#'                    'metric_d'=relative_error_delta)
#'     return(output)
#'   }
#' 


