#' CvM.normal
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
CvM.normal = function(x){
  pars <- estimate.normal(x)
  xbar <- pars[1]
  s <- pars[2]
  z <- pnorm(x,mean=xbar,sd=s)
  CvM(z)
}
