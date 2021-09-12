#' EDF statistics W^2 for Gamma Distribution
#'
#' Compute Cramer-von Mises statistic W^2 for an iid sample, x, to test for the Gamma distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.gamma" by default.
#'
#' @param x random sample
#' @param parameter parameters of Gamma distribution
#'
#' @return CvM.gamma gives Cramer-von Mises statistic of a uniform sample.
#' @export
#'
#' @examples
#' x = rgamma(100,1,1)
#' CvM.gamma(x)
#' CvM.gamma(x,parameter=c(1,1))
CvM.gamma <- function(x,parameter=estimate.gamma(x)){
    z <- cdf.gamma(x,parameter)
    CvM(z)
}
