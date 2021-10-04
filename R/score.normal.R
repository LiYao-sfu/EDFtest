#' Score Function for Normal Distribution
#'
#' Compute the score function for the two parameter Normal distribution
#'
#' @param x random sample
#' @param theta parameter of Normal distribution
#'
#' @return score.normal returns an n by 2 matrix whose ith row refers to the ith data point.
#'
#' @examples
score.normal = function(x,theta){
  sig=theta[2]
  mu=theta[1]
  s.mean= (x-mu)/sig
  s.sd= s.mean^2/sig-length(x)/sig
  cbind(s.mean/sig,s.sd)
}
