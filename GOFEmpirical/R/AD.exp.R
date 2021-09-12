#' EDF statistics A^2 for Exponential Distribution
#'
#' Compute Anderson-Darling statistic A^2 for an iid sample, x, to test for the Exponential distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.exp" by default.
#'
#' @param x random sample
#' @param parameter parameter of exponential distribution
#'
#' @return AD.exp gives Anderson-Darling statistic of a uniform sample.
#' @export
#'
#' @examples
#' x= rexp(100)
#' AD.exp(x)
#' AD.exp(x,1)
AD.exp = function(x,parameter=estimate.exp(x)){
  z = cdf.exp(x,parameter)
  AD(z)
}
