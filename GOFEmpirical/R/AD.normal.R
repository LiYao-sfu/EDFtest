#' AD.normal
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
AD.normal = function(x,parameter=estimate.normal(x)){
  #pars <- estimate.normal(x)
  xbar <- parameter[1]
  s <- parameter[2]
  z <- pnorm(x,mean=xbar,sd=s)
  AD(z)
}
