#' EDF statistics A^2
#'
#' Computes the Anderson-Darling goodness-of-fit statistic for testing uniformity on the unit interval.
#'
#' @param z Assumes z is a vector with 0 < z[i] < 1 for each i
#' @return float
#' @export
#'
#' @examples
#'
AD <- function(z){
  n <- length(z)
  u <- sort(z)
  i <- 2*(1:n)-1
  -n -sum(i*(log(u)+log(1-rev(u))))/n
}
