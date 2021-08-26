#' EDF statistics W^2
#'
#' Given probability integral transforms, compute the Cramer-von Mises goodness-of-fit statistic
#' for testing uniformity on the unit interval. It is a component of doing other tests.
#'
#' @param z vector of numbers supposed to be 0 and 1
#' @return CvM gives Cramer-von Mises statistic of a uniform sample.
#' @export
#'
#' @examples
#' x = runif(10)
#' CvM(x)
CvM <- function(z){
    n <- length(z)
    u <- sort(z)
    i <- (2*(1:n)-1)/(2*n)
    sum((u-i)^2)+1/(12*n)
}
