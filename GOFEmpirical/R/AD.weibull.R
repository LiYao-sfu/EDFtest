#' Anderson-Darling Statistics for the Weibull Distribution
#'
#' Computes Anderson-Darling statistic A^2 for an iid sample, x, to test for the Weibull distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.weibull"
#'
#' @param x numerical vector
#' @param parameter
#'
#' @return float
#' @export
#'
#' @examples
#'
AD.weibull <- function(x,parameter=estimate.weibull(x)){
#    pars <- estimate.weibull(x)
    alpha <- parameter[1]
    beta <- parameter[2]
    z <- pweibull(x,shape=alpha,scale=beta)
    AD(z)
}
