#' EDF statistics W^2 for Exponential Distribution
#'
#' Compute Cramer-von Mises statistic W^2 for an iid sample, x, to test for the Exponential distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.exp" by default.
#'
#' @param x random sample
#' @param parameter parameter of Exponential distribution
#'
#' @return CvM.exp gives Cramer-von Mises statistic of a uniform sample.
#' @export
#'
#' @examples
#' x= rexp(100)
#' CvM.exp(x)
#' CvM.exp(x,parameter=2)
CvM.exp = function(x,parameter=estimate.exp(x)){
  z = cdf.exp(x,parameter)
  CvM(z)
}
