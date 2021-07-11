#' cdf.gamma.regression
#'
#' @param y
#' @param x
#' @param theta.hat
#'
#' @return
#' @export
#'
#' @examples
cdf.gamma.regression = function(y,x,theta.hat){
  #
  #  Computes the probability integral transforms of the responses y
  #   in a gamma log-linear regression model: log(E(y)) = x beta
  #
  #  The mean of a Gamma distribution is shape * scale
  #
  coeff.hat=theta.hat[-1]
  shape.hat=theta.hat[1]
  linearpredictor = x%*%coeff.hat
  yhat=exp(linearpredictor)
  pgamma(y,shape=shape.hat,scale=yhat/shape.hat)
}
