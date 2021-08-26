#' EDF statistics A^2 for Normal Distribution
#'
#' Compute Anderson-Darling statistic A^2 for an iid sample, x, to test for the Normal distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.normal" by default.
#'
#' @param x random sample
#' @param parameter parameter of Normal distribution
#'
#' @return AD.normal gives Anderson-Darling statistic of a uniform sample.
#' @export
#'
#' @examples
#' x= rnorm(1000)
#' AD.normal(x)
#' AD.normal(x,c(0,1))
AD.normal = function(x,parameter=estimate.normal(x)){
  xbar <- parameter[1]
  s <- parameter[2]
  z <- pnorm(x,mean=xbar,sd=s)
  AD(z)
}
