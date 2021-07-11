#' score.laplace
#'
#' @param x
#' @param theta
#'
#' @return
#' @export
#'
#' @examples
score.laplace = function(x,theta){
  #
  # This is the score function for the two parameter Normal distribution
  #  It returns an n by 2 matrix whose ith row refers to the ith
  #  data point.
  #
  sig=theta[2]
  mu=theta[1]
  s.mean= signum(x-mu)/sig
  s.sd= s.mean^2/sig-length(x)/sig
  cbind(s.mean,s.sd)
}
