#' Score Function for Gamma Distribution
#'
#' Compute the score function for the two parameter Gamma distribution
#'
#' @param x random sample
#' @param theta parameter of Gamma distribution
#'
#' @return score.gamma returns an n by 2 matrix whose ith row refers to the ith data point.
score.gamma=function(x,theta){
  scale=theta[2]
  shape=theta[1]
  s.shape= log(x/scale)-digamma(shape)
  s.scale= x/scale^2 -shape/scale
  cbind(s.shape,s.scale)
}
