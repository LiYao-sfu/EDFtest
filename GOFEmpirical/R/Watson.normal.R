#' Watson.normal
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
Watson.normal = function(x){
  pars <- estimate.normal(x)
  xbar <- pars[1]
  s <- pars[2]
  z <- pnorm(x,mean=xbar,sd=s)
  Watson(z)
}
