#' Anderson-Darling Statistics for the Gamma Distribution
#'
#' Computes Anderson-Darling statistic A^2 for an iid sample, x, to test for the Gamma distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.gamma"
#'
#' @param x numerical vector
#' @param parameter
#'
#' @return float
#' @export
#'
#' @examples
#'
AD.gamma <- function(x,parameter=estimate.gamma(x)){
#    pars <- estimate.gamma(x)
    alpha <- parameter[1]
    beta <- parameter[2]
    z <- pgamma(x,shape=alpha,scale=beta)
    AD(z)
}
