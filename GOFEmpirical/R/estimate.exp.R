#' MLE for Exponential Distribution
#'
#' Computes the mle of the scale parameter for a sample x from the Exponential distribution.
#'
#' @param x random sample from the Exponential distirbution
#'
#' @return
#' @export
#'
#' @examples
#'
estimate.exp = function(x){
  mean(x)
}
