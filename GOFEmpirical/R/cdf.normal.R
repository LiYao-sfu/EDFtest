#' cdf.normal
#'
#' @param x
#' @param theta
#'
#' @return
#' @export
#'
#' @examples
cdf.normal = function(x,theta){
  #   shape=theta[1]
  #   scale=theta[2]
  pnorm(x,mean=theta[1],sd=theta[2])
}
