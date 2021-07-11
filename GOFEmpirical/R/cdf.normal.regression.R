#' cdf.normal.regression
#'
#' @param y
#' @param x
#' @param theta.hat
#'
#' @return
#' @export
#'
#' @examples
cdf.normal.regression = function(y,x,theta.hat){
  #
  #  Computes the probability integral transforms of the responses y
  #   in a normal linear regression model
  #
  coeff.hat=theta.hat[-1]
  sd.hat=theta.hat[1]
  y.hat = x%*%coeff.hat
  pnorm(y,mean=yhat,sd=sd.hat)
}
