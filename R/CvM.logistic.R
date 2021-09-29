#' EDF statistics W^2 for Logistic Distribution
#'
#' Compute Cramer-von Mises statistic W^2 for an iid sample, x, to test for the Logistic distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.logistic" by default.
#'
#' @param x random sample
#' @param parameter parameter of Logistic distribution
#'
#' @return CvM.logistic gives Cramer-von Mises statistic of a uniform sample.
#' @export
#'
#' @examples
#' x= rlogis(1000)
#' CvM.logistic(x)
#' CvM.logistic(x,c(0,1))
CvM.logistic <- function(x,parameter=estimate.logistic(x)){
    z <- cdf.logistic(x,parameter)
    CvM(z)
}
