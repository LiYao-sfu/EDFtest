#' Probability Integral Transforms for Normal Regression
#'
#' Computes the probability integral transforms of the responses y
#' in a normal linear regression model
#'
#' @param y response variable
#' @param x explanatory variables
#' @param theta.hat estimates of coefficients
#'
#' @return cdf.normal.regression gives probability integral transforms for y in a ordinary linear regression model
#' @export
#'
#' @examples
cdf.normal.regression = function(y,x,theta.hat){
  coeff.hat=theta.hat[-1]
  sd.hat=theta.hat[1]
  y.hat = x%*%coeff.hat
  pnorm(y,mean=y.hat,sd=sd.hat)
}
