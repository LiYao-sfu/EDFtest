#' EDF statistics W^2 for Laplace Distribution
#'
#' Compute Cramer-von Mises statistic W^2 for an iid sample, x, to test for the Laplace distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.laplace" by default.
#'
#' @param x random sample
#' @param parameter parameter of Laplace distribution
#'
#' @return CvM.laplace gives Cramer-von Mises statistic of a uniform sample.
#' @export
#'
#' @examples
#' library(L1pack)
#' x= rlaplace(1000,0,1)
#' estimate.laplace(x)
#' CvM.laplace(x)
#' CvM.laplace(x,c(0,1))
CvM.laplace = function(x,parameter=estimate.laplace(x)){
  z = cdf.laplace((x-parameter[1])/parameter[2])
  CvM(z)
}
