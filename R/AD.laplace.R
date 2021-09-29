#' EDF statistics A^2 for Laplace Distribution
#'
#' Compute Anderson-Darling statistic A^2 for an iid sample, x, to test for the Laplace distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.laplace" by default.
#'
#' @param x random sample
#' @param parameter parameter of Laplace distribution
#'
#' @return AD.laplace gives Anderson-Darling statistic of a uniform sample.
#' @export
#'
#' @examples
#' library(rmutil)
#' x= rlaplace(100)
#' estimate.laplace(x)
#' AD.laplace(x)
#' AD.laplace(x,c(0,1))
AD.laplace = function(x,parameter=estimate.laplace(x)){
  z = cdf.laplace(x,parameter)
  AD(z)
}
