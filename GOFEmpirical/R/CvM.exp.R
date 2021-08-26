#' EDF statistics W^2 for Exponential Distribution
#'
#' Compute Cramer-von Mises statistic W^2 for an iid sample, x, to test for the Gamma distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.gamma" by default.
#'
#' @param x random sample
#' @param parameter parameter of exponential distribution
#'
#' @return CvM.exp gives Cramer-von Mises statistic of a uniform sample.
#' @export
#'
#' @examples
#' x= rexp(100)
#' CvM.exp(x)
#' CvM.exp(x,parameter=2)
CvM.exp = function(x,parameter=estimate.exp(x)){
  theta = parameter[1]
  z = pexp(x/theta[1],rate=1)
  CvM(z)
}
