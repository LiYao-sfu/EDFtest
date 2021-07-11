#' score.weibull
#'
#' @param x
#' @param theta
#'
#' @return
#' @export
#'
#' @examples
score.weibull=function(x,theta){
  #
  # This is the score function for the two parameter Weibull distribution
  #
  scale=theta[2]
  shape=theta[1]
  r=x/scale
  lr=log(r)
  s.shape= 1/shape+lr-lr*r^shape
  s.scale= shape*(r^shape-1)/scale
  rbind(s.shape,s.scale)
}
