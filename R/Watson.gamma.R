#' EDF statistics U^2 for Gamma Distribution
#'
#' Compute Watson statistic U^2 for an iid sample, x, to test for the Gamma distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.gamma" by default.
#'
#' @param x random sample
#' @param parameter parameter of Gamma distribution
#'
#' @return Watson.gamma gives Watson statistic of a uniform sample
#' @export
#'
#' @examples
#' x= rgamma(1000,1)
#' estimate.gamma(x)
#' Watson.gamma(x)
Watson.gamma = function(x,parameter=estimate.gamma(x)){
  z <- cdf.gamma(x,parameter)
  Watson(z)
}
