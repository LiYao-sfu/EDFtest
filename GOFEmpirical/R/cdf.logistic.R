#' Probability Integral Transforms for Logistic Distribution
#'
#' Compute the probability integral transforms of the random sample x in a Logistic distribution
#'
#' @param x random sample
#' @param theta parameters of Logistic distribution
#'
#' @return cdf.logistic gives probability integral transforms for dataset following Logistic distribution
#' @export
#'
#' @examples
#' x= rlogis(100)
#' cdf.logistic(x,c(0,1))
cdf.logistic = function(x,theta){
  plogis(x,location=theta[1],scale=theta[2])
}
