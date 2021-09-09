#' Probability Integral Transforms for Exponential Distribution
#'
#' Compute the probability integral transforms of the random sample x in a Exponential distribution
#'
#' @param x random sample
#' @param theta parameter of Exponential distribution
#'
#' @return cdf.exp gives probability integral transforms for dataset following Weibull distribution
#' @export
#'
#' @examples
#' x= rexp(100,1)
#' cdf.exp(x,1)
cdf.exp = function(x,theta){
  pexp(x,rate=1/theta)
}
