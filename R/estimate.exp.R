#' MLE for Exponential Distribution
#'
#' Estimate scale parameter of the Exponential distribution by the method of maximum likelihood.
#'
#' @param x random sample
#' @param use.rate Logical, if TRUE use the rate instead of the scale
#'
#' @return estimated scale (or rate) parameter of the Exponential distribution.
#' @export
#'
#' @examples
#' x=rexp(10)
#' estimate.exp(x)
estimate.exp = function(x,use.rate=FALSE){
  if(use.rate){
    return(1/mean(x))
  }
  mean(x)
}
