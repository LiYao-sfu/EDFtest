#' EDF statistics U^2
#'
#' Given probability integral transforms, compute the Watson goodness-of-fit statistic
#' for testing uniformity on the unit interval. It is also a component of doing other tests.
#'
#' @param z vector of numbers supposed to be 0 and 1.
#'
#' @return Watson gives Watson statistic of a uniform sample.
#' @export
#'
#' @examples
#' x = runif(10)
#' Watson(x)
Watson <- function(z){
    n <- length(z)
    u <- sort(z)
    i <- (2*(1:n)-1)/(2*n)
    w=sum((u-i)^2)+1/(12*n)
    corr=n*(mean(u)-0.5)^2 #correction
    w-corr
}
