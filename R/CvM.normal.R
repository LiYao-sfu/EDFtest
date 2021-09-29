#' EDF statistics W^2 for Normal Distribution
#'
#' Compute Cramer-von Mises statistic W^2 for an iid sample, x, to test for the Normal distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.normal" by default.
#'
#' @param x random sample
#' @param parameter parameter of Normal distribution
#'
#' @return CvM.normal gives Cramer-von Mises statistic of a uniform sample.
#' @export
#'
#' @examples
#' x= rnorm(1000)
#' CvM.normal(x)
#' CvM.normal(x,c(0,1))
CvM.normal = function(x,parameter=estimate.normal(x)){
  z <- cdf.normal(x,parameter)
  CvM(z)
}
