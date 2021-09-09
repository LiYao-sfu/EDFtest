#' Probability Integral Transforms for Weibull Distribution
#'
#' Compute the probability integral transforms of the random sample x in a Weibull distribution
#'
#' @param x random sample
#' @param theta parameters of Weibull distribution
#'
#' @return cdf.weibull gives probability integral transforms for dataset following Weibull distribution
#' @export
#'
#' @examples
#' x= rweibull(100,1)
#' cdf.weibull(x,c(1,1))
cdf.weibull = function(x,theta){
  pweibull(x,shape=theta[1],scale=theta[2])
}
