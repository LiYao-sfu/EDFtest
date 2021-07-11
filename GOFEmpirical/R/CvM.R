#' EDF statistics W^2
#'
#' Computes the Cramer-von Mises goodness-of-fit statistic for testing uniformity on the unit interval.
#'
#' @param z Assumes z is a vector with 0 < z[i] < 1 for each i
#' @return float
#' @export
#'
#' @examples
#'
CvM <- function(z){
    n <- length(z)
    u <- sort(z)
    i <- (2*(1:n)-1)/(2*n)
    sum((u-i)^2)+1/(12*n)
  }
