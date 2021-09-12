#' EDF statistics U^2 for Weibull Distribution
#'
#' Compute Watson statistic U^2 for an iid sample, x, to test for the Weibull distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.weibull" by default.
#'
#' @param x random sample
#' @param parameter parameter of Weibull distribution
#'
#' @return Watson.weibull gives Watson statistic of a uniform sample
#' @export
#'
#' @examples
#' x= rweibull(1000,1)
#' estimate.weibull(x)
#' Watson.weibull(x)
Watson.weibull = function(x,parameter=estimate.weibull(x)){
  z <- cdf.weibull(x,parameter)
  Watson(z)
}
