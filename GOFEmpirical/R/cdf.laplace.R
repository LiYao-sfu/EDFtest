#' cdf.laplace
#'
#' @param x
#' @param theta
#'
#' @return
#' @export
#'
#' @examples
cdf.laplace = function(x,theta=c(0,1)){
  #   shape=theta[1]
  #   scale=theta[2]
  m = theta[1]
  s = theta[2]
  v = exp(x-m)/2
  v[x>m] = 1-exp(-(x-m)[x>m])/2
  v
}
