#' Cramer Von Mises Statistics for the Logistic Distribution
#'
#' Computes Cramer Von Mises statistic W^2 for an iid sample, x, to test for the Logistic distribution with parameters unknown.
#' Estimate parameters by default ML using "estimate.logistic"
#'
#' @param x
#' @param parameter
#'
#' @return
#' @export
#'
#' @examples
#'
CvM.logistic <- function(x,parameter=estimate.logistic(x)){
#    pars <- estimate.logistic(x)
    alpha <- parameter[1]
    beta <- parameter[2]
    z <- plogis(x,location=alpha,scale=beta)
    CvM(z)
}
