#' Watson Statistics
#'
#' @param z
#'
#' @return
#' @export
#'
#' @examples
Watson <- function(z){
    n <- length(z)
    u <- sort(z)
    i <- (2*(1:n)-1)/(2*n)
    w=sum((u-i)^2)+1/(12*n)
    corr=n*(mean(u)-0.5)^2 #correction
    w-corr
}
