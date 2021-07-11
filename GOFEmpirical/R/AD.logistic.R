#' Anderson-Darling Statistics for the Logistic Distribution
#'
#' Computes Anderson-Darling statistic A^2 for an iid sample, x, to test for the Logistic distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.logistic"
#'
#' @param x numerical vector
#' @param parameter
#'
#' @return float
#' @export
#'
#' @examples
#'
AD.logistic <- function(x,parameter=estimate.logistic(x)){
#    pars <- estimate.logistic(x)
    alpha <- parameter[1]
    beta <- parameter[2]
    z <- plogis(x,location=alpha,scale=beta)
    AD(z)
}
