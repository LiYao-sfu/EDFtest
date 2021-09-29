#' EDF statistics A^2
#'
#' Given probability integral transforms, compute the Anderson-Darling goodness-of-fit statistic
#' for testing uniformity on the unit interval.
#'
#' @param z vector of numbers supposed to be 0 and 1
#' @return AD gives Anderson-Darling statistic of a uniform sample.
AD <- function(z){
  n <- length(z)
  u <- sort(z)
  i <- 2*(1:n)-1
  -n -sum(i*(log(u)+log(1-rev(u))))/n
}
