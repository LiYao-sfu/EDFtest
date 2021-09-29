#' MLE for Laplace (double exponential) Distribution
#'
#' Estimate location and scale parameters of the Laplace distribution by the method of maximum likelihood.
#'
#' @param x random sample
#' @param use.sd Logical, if TRUE use the sd instead of the scale
#'
#' @return estimated location and scale (or sd) parameters of the Laplace distribution.
#' @export
#'
#' @examples
#' library(rmutil)
#' library(L1pack)
#' x= rmutil::rlaplace(1000,0,1)
#' estimate.laplace(x)
#' x= L1pack::rlaplace(1000,0,1)
#' estimate.laplace(x,use.sd=TRUE)
estimate.laplace = function(x,use.sd=FALSE){
  med = median(x)
  MAD = mean(abs(x-med))
  if(use.sd){
    return(c(med,sqrt(2)*MAD))
  }
  c(med,MAD)
}
