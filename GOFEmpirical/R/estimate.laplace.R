#' MLE for Laplace (double exponential) Distribution
#'
#' Computes the mle of the location, mu, and scale, b, for a sample x from the Laplace distribution.
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
#'
estimate.laplace = function(x){
  med = median(x)
  MAD = mean(abs(x-med))
  c(med,MAD)
}
