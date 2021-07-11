#' cdf.gamma
#'
#' @param x
#' @param theta
#'
#' @return
#' @export
#'
#' @examples
cdf.gamma = function(x,theta){
  #   shape=theta[1]
  #   scale=theta[2]
  pgamma(x,shape=theta[1],scale=theta[2])
}
