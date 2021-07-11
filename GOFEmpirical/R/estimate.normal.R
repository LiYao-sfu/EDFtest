#' MLE for Normal Distribution
#'
#' Computes the mean and variance of a sample x from the normal distribution.
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
estimate.normal = function(x){
  c(mean(x),sd(x))

}
