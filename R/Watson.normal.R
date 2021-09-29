#' EDF statistics U^2 for Normal Distribution
#'
#' Compute Watson statistic U^2 for an iid sample, x, to test for the Normal distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.normal" by default.
#'
#' @param x random sample
#' @param parameter parameter of Normal distribution
#'
#' @return Watson.normal gives Watson statistic of a uniform sample.
#' @export
#'
#' @examples
#' x= rnorm(1000)
#' Watson.normal(x)
#' Watson.normal(x,c(0,1))
Watson.normal = function(x,parameter=estimate.normal(x)){
  z <- cdf.normal(x, theta=parameter)
  Watson(z)
}
