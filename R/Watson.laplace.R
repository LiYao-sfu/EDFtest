#' EDF statistics U^2 for Laplace Distribution
#'
#' Compute Watson statistic U^2 for an iid sample, x, to test for the Laplace distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.laplace" by default.
#'
#' @param x random sample
#' @param parameter parameter of Laplace distribution
#'
#' @return Watson.laplace gives Watson statistic of a uniform sample
#' @export
#'
#' @examples
#' library(rmutil)
#' x= rlaplace(1000,0,1)
#' estimate.laplace(x)
#' Watson.laplace(x,c(0,1))
Watson.laplace = function(x,parameter=estimate.laplace(x)){
  z <- cdf.laplace(x,parameter)
  Watson(z)
}
