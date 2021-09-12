#' EDF statistics U^2 for Exponential Distribution
#'
#' Compute Watson statistic U^2 for an iid sample, x, to test for the Exponential distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.exp" by default.
#'
#' @param x random sample
#' @param parameter parameter of Exponential distribution
#'
#' @return Watson.exp gives Watson statistic of a uniform sample
#' @export
#'
#' @examples
#' x= rexp(1000,2)
#' estimate.exp(x)
#' Watson.exp(x)
Watson.exp = function(x,parameter=estimate.exp(x)){
  z <- cdf.exp(x,parameter)
  Watson(z)
}
