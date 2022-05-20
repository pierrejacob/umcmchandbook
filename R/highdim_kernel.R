### TAKEN FROM https://github.com/niloyb/CoupledHalfT
### and slightly modified

## Matrix calculation helper functions ##
#' X_eta_tX
#' @description Calculates X * Diag(1/eta) * t(X_transpose)
#' @param eta vector of length p
#' @param X matrix of length n by p
#' @param X_transpose Pre-calculated transpose of X
#' @return X Diag(1/eta) t(X_transpose)
#' @export
X_eta_tX <- function(eta, X, X_transpose){
  # NOTE: (1) For MacOS with veclib BLAS, crossprod is fast via multi-threading
  # return(crossprod(X_transpose*c(1/eta)^0.5))
  return(fcprd(X_transpose * c(1/eta)^0.5))
}
#' M_matrix
#' @description Calculates I + X * Diag(1/eta) * t(X_transpose)
#' @param xi positive scalar
#' @param eta vector of length p
#' @param X_eta_tX_matrix X * Diag(1/eta) * t(X_transpose), matrix n by n
#' @param n positive integer
#' @return I + X Diag(1/eta) t(X_transpose)
#' @export
M_matrix <- function(xi, eta, X_eta_tX_matrix, n){
  if (length(eta) == 0) return(diag(n))
  return(diag(n) + (xi^-1)*X_eta_tX_matrix)
}

## xi update given eta ##
# Unnormalized posterior pdf of log(xi)
#' log_ratio
#' @description Unnormalized posterior pdf of log(xi)
#' @param xi positive scalar
#' @param eta length p vector
#' @param X_eta_tX_matrix X * Diag(1/eta) * t(X_transpose), matrix n by n
#' @param y length n vector
#' @param a0 positive scalar
#' @param b0 positive scalar
#' @return Unnormalized posterior pdf of log(xi)
#' @export
log_ratio <- function(xi, eta, X_eta_tX_matrix, y, a0, b0){
  n <- length(y)
  M <- M_matrix(xi,eta,X_eta_tX_matrix,n)
  chol_M <- chol(M)
  log_det_M <- 2*sum(log(diag(chol_M)))
  M_inverse <- chol2inv(chol_M)
  ssr <- b0 + t(y)%*%((M_inverse)%*%y)
  log_likelihood <- -0.5*log_det_M -0.5*(n+a0)*log(ssr)
  log_prob <- -log(sqrt(xi)*(1+xi))
  return(list('log_likelihood' = log_likelihood+log_prob,
              'ssr' = ssr, 'M_matrix_inverse' = M_inverse))
}

## Metropolis-Rosenbluth-Teller-Hastings update of xi given eta
#' xi_update
#' @description Metropolis-Hastings update of xi
#' @param xi positive scalar
#' @param eta length p vector
#' @param X n by p matrix
#' @param X_transpose Pre-calculated transpose of X
#' @param y length n vector
#' @param a0 positive scalar
#' @param b0 positive scalar
#' @param std_MH standard deviation of log-normal MH proposal
#' @return Metropolis-Hastings update of xi given eta
#' @export
xi_update <- function(xi, eta, X, X_transpose, y, a0, b0, std_MH){
  # n <- length(y)
  proposed_xi <- exp(rnorm(1, mean = log(xi), sd = std_MH))
  X_eta_tX_matrix <- X_eta_tX(eta, X, X_transpose)
  log_ratio_current_and_ssr <- log_ratio(xi, eta,
                                         X_eta_tX_matrix, y, a0, b0)
  log_ratio_proposed_and_ssr <- log_ratio(proposed_xi, eta,
                                          X_eta_tX_matrix, y, a0, b0)
  # MRTH accept-reject step
  log_accept_prob <- (log_ratio_proposed_and_ssr$log_likelihood - log_ratio_current_and_ssr$log_likelihood) +
    (log(proposed_xi) - log(xi))
  if (log(runif(1))<log_accept_prob){
    xi <- proposed_xi
    log_ratio_current_and_ssr <- log_ratio_proposed_and_ssr
  }
  return(list('xi'= xi, 'ssr' = log_ratio_current_and_ssr$ssr,
              'M_matrix' = log_ratio_current_and_ssr$M_matrix_inverse))
}

## sigma^2 update step given eta and xi
#' sigma2_update
#' @description sigma^2 update step given eta and xi
#' @param xi positive scalar
#' @param eta length p vector
#' @param n positive integer
#' @param ssr positive scalar
#' @param a0 positive scalar
#' @param b0 positive scalar
#' @return sigma^2 update step given eta and xi
#' @export
sigma2_update <- function(xi, eta, n, ssr, a0, b0){
  # if(!is.null(sigma2_fixed_value)) return(sigma2_fixed_value)
  return(1/(rgamma(1, shape = (n+a0)/2, rate = (ssr)/2)))
}

## beta update given xi, eta, sigma2
#' beta_update
#' @description beta given xi, eta, sigma2 update using algo of Bhattacharya
#' @param xi positive scalar
#' @param sigma2 positive scalar
#' @param eta length p vector
#' @param X n by p matrix
#' @param X_transpose Pre-calculated transpose of X
#' @param y length n vector
#' @param M_matrix_inverse n by n matrix
#' @return beta update given xi, eta, sigma2
#' @export
beta_update <- function(xi, sigma2, eta, X, X_transpose, y, M_matrix_inverse){
  p <- length(eta)
  n <- length(y)
  u = rnorm(p, 0, 1)
  u = (sqrt(xi*eta)^-1)*u
  v = X%*%u + c(rnorm(n,0,1))
  v_star <- M_matrix_inverse%*%(y/sqrt(sigma2) - v)
  U <- (xi^-1)* ( ((eta)^(-1))*(X_transpose) )
  u <- u + U%*%v_star
  beta <- sqrt(sigma2)*(u)
  return(beta)
}


## Full blocked Gibbs samplers ##
#' half_t_kernel
#' @description Full blocked Gibbs kernel for half-t priors
#' @param X n by p matrix
#' @param X_transpose Pre-calculated transpose of X
#' @param y length n vector
#' @param a0 positive scalar
#' @param b0 positive scalar
#' @param std_MH standard deviation of log-normal MH proposal
#' @param xi_current current xi value (positive scalar)
#' @param sigma2_current current sigma2 value (positive scalar)
#' @param beta_current current beta value (vector of length p)
#' @param eta_current current eta value (vector of length p)
#' @param nrepeats_eta number of slice sampling steps
#' @param t_dist_df degree of freedom v for Half-t(v). Take v >=1.
#' @return (beta, eta, sigma2, xi) sampled from the half-t prior blocked Gibbs kernel
#' @export
half_t_kernel <-
  function(X, X_transpose, y, a0=1, b0=1, std_MH=0.8,
           xi_current, sigma2_current,
           beta_current, eta_current,
           nrepeats_eta = 1, 
           t_dist_df){
  n <- dim(X)[1]
  p <- dim(X)[2]
  # eta update
  eta_new <- eta_update(xi_current, sigma2_current, beta_current, eta_current, t_dist_df = t_dist_df, nrepeats = nrepeats_eta)
  if (!is.null(dim(eta_new))) eta_new <- as.numeric(eta_new)
  # xi update
  xi_new <- xi_update(xi_current, eta_new, X, X_transpose, y, a0, b0, std_MH)
  # sigma2 update
  sigma2_new <- sigma2_update(xi_new$xi, eta_new, n, xi_new$ssr, a0, b0)
  # beta update
  beta_new <- beta_update(xi_new$xi, sigma2_new, eta_new, X, X_transpose, y,
                          xi_new$M_matrix)
  output <- list( 'beta'=beta_new, 'eta'=eta_new,
                  'sigma2'=sigma2_new, 'xi'=xi_new$xi)
  return(output)
}




