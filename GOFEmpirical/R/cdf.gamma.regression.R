#' Probability Integral Transforms for Gamma Regression
#'
#' Computes the probability integral transforms of the responses y
#' in a gamma log-linear regression model: log(E(y)) = x beta
#'
#' @param y response variable
#' @param x explanatory variables
#' @param theta.hat estimates of coefficients
#'
#' @return cdf.gamma.regression gives probability integral transforms for y in a gamma log-linear regression model
#' @export
#'
#' @examples
cdf.gamma.regression = function(y,x,theta.hat){
  coeff.hat=theta.hat[-1]
  shape.hat=theta.hat[1]
  linearpredictor = x%*%coeff.hat
  yhat=exp(linearpredictor)
  pgamma(y,shape=shape.hat,scale=yhat/shape.hat)
}
