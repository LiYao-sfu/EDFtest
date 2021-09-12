#' EDF statistics U^2 for Logistic Distribution
#'
#' Compute Watson statistic U^2 for an iid sample, x, to test for the Logistic distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.logistic" by default.
#'
#' @param x random sample
#' @param parameter parameter of Logistic distribution
#'
#' @return Watson.logistic gives Watson statistic of a uniform sample
#' @export
#'
#' @examples
#' x= rlogis(1000)
#' estimate.logistic(x)
#' Watson.logistic(x)
Watson.logistic = function(x,parameter=estimate.logistic(x)){
  z <- cdf.logistic(x,parameter)
  Watson(z)
}
