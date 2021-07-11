#' Cramer Von Mises Statistics for the Gamma Distribution
#'
#' Computes Cramer Von Mises statistic W^2 for an iid sample, x, to test for the Gamma distribution with parameters unknown.
#' Estimate parameters by default ML using "estimate.gamma"
#'
#' @param x
#' @param parameter
#'
#' @return
#' @export
#'
#' @examples
#'
CvM.gamma <- function(x,parameter=estimate.gamma(x)){
#    pars <- estimate.gamma(x)
    alpha <- parameter[1]
    beta <- parameter[2]
    z <- pgamma(x,shape=alpha,scale=beta)
    CvM(z)
}
