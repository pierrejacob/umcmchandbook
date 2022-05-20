### unbiased point-wise estimate of Fishy functions
## via the representation \tilde{h}(x) = sum_{n \geq 0} {P^n h (x) - P^n h (x_0)}
## obtained by sum_{n \geq 0} {h(X_n) - h(Y_n)} with X_n started at x and Y_n started at x_0
## state_x and state_x_0 should be list with entries 'position' and whatever else
## is needed by coupled_kernel
#'@rdname sample_unbiasedfishy
#'@title Sample unbiased estimator of evaluation of fishy function
#'@description Sample two Markov chains until they meet. One 
#'starts from x and the other from y. The output variable is an unbiased estimator
#' of a fishy function \eqn{x \mapsto \sum_{t\geq 0} \{ P^t h(x) - P^t h (y) \}}.
#' 
#'@param coupled_kernel A list taking two states and returning two states, performing one step of a coupled Markov kernel;
#'it also returns a boolean "identical" indicating whether the two states are identical.
#'@param h test function with which we want to solve the Poisson equation
#'@param state_x state of Markov chain at position x
#'@param state_x_0 state of Markov chain at position y
#'@param max_iterations A maximum number of iterations, at which to interrup the while loop; Inf by default
#'@return A list with
#'\itemize{
#'
#'\item meetingtime: the meeting time; equal to Inf if while loop was interrupted
#'
#'\item estimator: value of the unbiased estimator of fishy function at x
#'}
#'@export
sample_unbiasedfishy <- function(coupled_kernel, h, state_x, state_x_0, max_iterations = Inf){
  meetingtime <- Inf
  time <- 0
  estimator_ <- h(state_x$position) - h(state_x_0$position)
  while (is.infinite(meetingtime) && (time < max_iterations)){
    time <- time + 1
    # use coupled kernel
    coupledstates <- coupled_kernel(state_x, state_x_0)
    state_x <- coupledstates$state1
    state_x_0 <- coupledstates$state2
    estimator_ <- estimator_ + h(state_x$position) - h(state_x_0$position)
    # check if meeting happens
    if (coupledstates$identical) meetingtime <- time
  }
  return(list(estimator = estimator_, meetingtime = meetingtime))
}

