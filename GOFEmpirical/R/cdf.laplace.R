#' Probability Integral Transforms for Laplace Distribution
#'
#' Compute the probability integral transforms of the random sample x in a Laplace distribution
#'
#' @param x random sample
#' @param theta parameters of Laplace distribution
#'
#' @return cdf.laplace gives probability integral transforms for dataset following Laplace distribution
#' @export
#'
#' @examples
#' library(L1pack)
#' x= rlaplace(100)
#' cdf.laplace(x,c(0,1))
cdf.laplace = function(x,theta){
  m = theta[1]
  s = theta[2]
  v = exp(x-m)/2
  v[x>m] = 1-exp(-(x-m)[x>m])/2
  v
}
