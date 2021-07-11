#' Cramer Von Mises Statistics for the Weibull Distribution
#'
#' Computes Cramer Von Mises statistic W^2 for an iid sample, x, to test for the Weibull distribution with parameters unknown.
#' Estimate parameters by default ML using "estimate.weibull"
#'
#' @param x
#' @param parameter
#'
#' @return
#' @export
#'
#' @examples
#'
CvM.weibull <- function(x,parameter=estimate.weibull(x)){
#    pars <- estimate.weibull(x)
    alpha <- parameter[1]
    beta <- parameter[2]
    z <- pweibull(x,shape=alpha,scale=beta)
    CvM(z)
}
