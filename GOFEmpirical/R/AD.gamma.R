#' EDF statistics A^2 for Gamma Distribution
#'
#' Compute Anderson-Darling statistic A^2 for an iid sample, x, to test for the Gamma distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.gamma" by default.
#'
#' @param x random sample
#' @param parameter parameters of Gamma distribution
#'
#' @return AD.gamma gives Anderson-Darling statistic of a uniform sample.
#' @export
#'
#' @examples
#' x = rgamma(100,1,2)
#' AD.gamma(x)
#' AD.gamma(x,c(1,1/2))
AD.gamma <- function(x,parameter=estimate.gamma(x)){
    alpha <- parameter[1]
    beta <- parameter[2]
    z <- pgamma(x,shape=alpha,scale=beta)
    AD(z)
}
