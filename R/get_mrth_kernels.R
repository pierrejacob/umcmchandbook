#'@rdname get_mrth_kernel
#'@title Get random walk Metropolis-Rosenbluth-Teller-Hastings kernels for 1d target
#'@description This function takes a target (specified through its log-pdf)
#' and a variance for a Normal random walk proposal, and returns a list containing the keys
#' \code{single_kernel}, \code{coupled_kernel} corresponding to marginal
#' and coupled MH kernels.
#'
#' The coupling is done by reflection-maximal coupling of the proposals,
#' and common uniform variable for the accept/reject step. For reflection-maximal
#' couplings, see \code{\link{rnorm_reflectionmax}}.
#'
#' The returned kernels can then be used in the functions \code{\link{sample_meetingtime}} or
#' \code{\link{sample_coupled_chains}} or \code{\link{sample_unbiasedestimator}}.
#'
#'@param target function taking a vector as input and returning target log-density evaluation
#'@param Sigma_proposal variance of the Normal random walk proposal
#'@return A list containing the keys
#' \code{single_kernel}, \code{coupled_kernel}.
#'@export
get_mrth_kernels <- function(target, Sigma_proposal){
  dimension <- 1
  Sigma_proposal_chol <- sqrt(Sigma_proposal)
  zeromean <- 0
  # single kernel
  single_kernel <- function(state){
    position <- state$position
    current_pdf <- state$current_pdf
    proposal_value <- rnorm(1, position, Sigma_proposal_chol)
    proposal_pdf <- target(proposal_value)
    accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
    if (accept){
      return(list(position = proposal_value, current_pdf = proposal_pdf, accept = accept))
    } else {
      return(list(position = position, current_pdf = current_pdf, accept = accept))
    }
  }
  # coupled kernel
  coupled_kernel <- function(state1, state2){
    position1 <- state1$position; current_pdf1 <- state1$current_pdf
    position2 <- state2$position; current_pdf2 <- state2$current_pdf
    proposal_value <- rnorm_reflectionmax(position1, position2, Sigma_proposal_chol)
    proposal1 <- proposal_value$xy[1]; proposal_pdf1 <- target(proposal1)
    if (proposal_value$identical){
      proposal2 <- proposal1; proposal_pdf2 <- proposal_pdf1
    } else {
      proposal2 <- proposal_value$xy[2]; proposal_pdf2 <- target(proposal2)
    }
    logu <- log(runif(1))
    accept1 <- FALSE; accept2 <- FALSE
    if (is.finite(proposal_pdf1)){
      accept1 <- (logu < (proposal_pdf1 - current_pdf1))
    }
    if (is.finite(proposal_pdf2)){
      accept2 <- (logu < (proposal_pdf2 - current_pdf2))
    }
    if (accept1){
      position1 <- proposal1
      current_pdf1 <- proposal_pdf1
    }
    if (accept2){
      position2 <- proposal2
      current_pdf2 <- proposal_pdf2
    }
    identical_ <- proposal_value$identical && accept1 && accept2
    return(list(state1 = list(position = position1, current_pdf = current_pdf1),
                state2 = list(position = position2, current_pdf = current_pdf2),
                identical = identical_))
  }
  return(list(single_kernel = single_kernel, coupled_kernel = coupled_kernel))
}
