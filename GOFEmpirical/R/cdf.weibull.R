#' cdf.weibull
#'
#' @param x
#' @param theta
#'
#' @return
#' @export
#'
#' @examples
cdf.weibull = function(x,theta){
  #   shape=theta[1]
  #   scale=theta[2]
  pweibull(x,shape=theta[1],scale=theta[2])
}
