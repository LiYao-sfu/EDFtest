#' Probability Integral Transforms for Gamma Distribution
#'
#' Compute the probability integral transforms of the random sample x in a Gamma distribution
#'
#' @param x random sample
#' @param theta parameters of Gamma distribution
#'
#' @return cdf.gamma gives probability integral transforms for dataset following Gamma distribution
#' @export
#'
#' @examples
#' x = rgamma(100,1,2)
#' cdf.gamma(x,c(1,1/2))
cdf.gamma = function(x,theta){
  pgamma(x,shape=theta[1],scale=theta[2])
}
