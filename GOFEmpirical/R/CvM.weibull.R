#' EDF statistics W^2 for Weibull Distribution
#'
#' Compute Cramer-von Mises statistic W^2 for an iid sample, x, to test for the Weibull distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.weibull" by default.
#'
#' @param x random sample
#' @param parameter parameter of Weibull distribution
#'
#' @return CvM.weibull gives Cramer-von Mises statistic of a uniform sample.
#' @export
#'
#' @examples
#' x= rweibull(1000,1)
#' CvM.weibull(x)
#' CvM.weibull(x,c(1,1))
CvM.weibull <- function(x,parameter=estimate.weibull(x)){
    alpha <- parameter[1]
    beta <- parameter[2]
    z <- pweibull(x,shape=alpha,scale=beta)
    CvM(z)
}
