#' Probability Integral Transforms for Normal Distribution
#'
#' Compute the probability integral transforms of the random sample x in a Normal distribution
#'
#' @param x random sample
#' @param theta parameters of Normal distribution
#'
#' @return cdf.logistic gives probability integral transforms for dataset following Normal distribution
#' @export
#'
#' @examples
#' x= rnorm(100)
#' cdf.normal(x,c(0,1))
cdf.normal = function(x,theta){
  pnorm(x,mean=theta[1],sd=theta[2])
}
