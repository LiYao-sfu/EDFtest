#' Score Function for Weibull Distribution
#'
#' Compute the score function for the two parameter Weibull distribution
#'
#' @param x random sample
#' @param theta parameter of Weibull distribution
#'
#' @return score.weibull returns an n by 2 matrix whose ith row refers to the ith data point.
#'
#' @examples
score.weibull=function(x,theta){
  scale=theta[2]
  shape=theta[1]
  r=x/scale
  lr=log(r)
  s.shape= 1/shape+lr-lr*r^shape
  s.scale= shape*(r^shape-1)/scale
  cbind(s.shape,s.scale)
}
