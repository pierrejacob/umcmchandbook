#'@export
tv_upper_bound <- function(meetingtimes, lag, t){
  return(mean(pmax(0, ceiling((meetingtimes-lag-t) / lag))))
}
