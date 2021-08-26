#' EDF statistics A^2 for Weibull Distribution
#'
#' Compute Anderson-Darling statistic A^2 for an iid sample, x, to test for the Weibull distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.weibull" by default.
#'
#' @param x random sample
#' @param parameter parameter of Weibull distribution
#'
#' @return AD.weibull gives Anderson-Darling statistic of a uniform sample.
#' @export
#'
#' @examples
#' x= rweibull(1000,1)
#' AD.weibull(x)
#' AD.weibull(x,c(1,1))
AD.weibull <- function(x,parameter=estimate.weibull(x)){
    alpha <- parameter[1]
    beta <- parameter[2]
    z <- pweibull(x,shape=alpha,scale=beta)
    AD(z)
}
