#' Score Function for Laplace Distribution
#'
#' Compute the score function for the two parameter Laplace distribution
#'
#' @param x random sample
#' @param theta parameter of Laplace distribution
#'
#' @return score.laplace returns an n by 2 matrix whose ith row refers to the ith data point.
#'
#' @examples
score.laplace = function(x,theta){
  sig=theta[2]
  mu=theta[1]
  signum<-function(x){
    y=x/abs(x)
    y[is.na(y)]=0
    return(y)
  }
  s.mean= signum(x-mu)/sig
  s.sd= s.mean^2/sig-length(x)/sig
  cbind(s.mean,s.sd)
}
