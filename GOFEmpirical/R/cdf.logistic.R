#' cdf.logistic
#'
#' @param x
#' @param theta
#'
#' @return
#' @export
#'
#' @examples
cdf.logistic = function(x,theta){
  #   shape=theta[1]
  #   scale=theta[2]
  plogis(x,location=theta[1],scale=theta[2])
}
