#' MLE for Normal Distribution
#'
#' Estimate location and scale parameters of the Normal distribution by the method of maximum likelihood.
#'
#' @param x random sample
#'
#' @return estimated location and scale parameters of the Normal distribution.
#' @export
#'
#' @examples
#' x=rnorm(10)
#' estimate.normal(x)
estimate.normal = function(x){
  c(mean(x),sd(x))

}
