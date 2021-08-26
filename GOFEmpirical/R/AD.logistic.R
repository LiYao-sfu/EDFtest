#' EDF statistics A^2 for Logistic Distribution
#'
#' Compute Anderson-Darling statistic A^2 for an iid sample, x, to test for the Logistic distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.logistic" by default.
#'
#' @param x random sample
#' @param parameter parameter of Logistic distribution
#'
#' @return AD.logistic gives Anderson-Darling statistic of a uniform sample.
#' @export
#'
#' @examples
#' x= rlogis(1000)
#' AD.logistic(x)
#' AD.logistic(x,c(0,1))
AD.logistic <- function(x,parameter=estimate.logistic(x)){
    alpha <- parameter[1]
    beta <- parameter[2]
    z <- plogis(x,location=alpha,scale=beta)
    AD(z)
}
