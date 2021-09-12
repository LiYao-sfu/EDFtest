#' MLE for Laplace (double exponential) Distribution
#'
#' Estimate location and scale parameters of the Laplace distribution by the method of maximum likelihood.
#'
#' @param x random sample
#'
#' @return estimated location and scale parameters of the Laplace distribution.
#' @export
#'
#' @examples
#' library(rmutil)
#' x= rlaplace(1000,0,1)
#' estimate.laplace(x)
estimate.laplace = function(x){
  med = median(x)
  MAD = mean(abs(x-med))
  c(med,MAD)
}
