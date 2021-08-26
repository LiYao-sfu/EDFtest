#' MLE for Exponential Distribution
#'
#' Estimate scale parameter of the Exponential distribution by the method of maximum likelihood.
#'
#' @param x random sample
#'
#' @return estimated scale parameter of the Exponential distribution.
#' @export
#'
#' @examples
#' x=rexp(10)
#' estimate.exp(x)
estimate.exp = function(x){
  mean(x)
}
