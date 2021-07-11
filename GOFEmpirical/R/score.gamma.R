#' score.gamma
#'
#' @param x
#' @param theta
#'
#' @return
#' @export
#'
#' @examples
score.gamma=function(x,theta){
  #
  # This is the score function for the two parameter Gamma distribution
  #
  scale=theta[2]
  shape=theta[1]
  s.shape= log(x/scale)-digamma(shape)
  s.scale= x/scale^2 -shape/scale
  rbind(s.shape,s.scale) #should be cbind
}
