#' Probability Integral Transforms for Laplace Regression
#'
#' Computes the probability integral transforms of the responses y
#' in a normal linear regression model
#'
#' @param y response variable
#' @param x explanatory variables
#' @param theta.hat estimates of coefficients
#'
#' @return cdf.laplace.regression gives probability integral transforms for y in a
#' @export
#'
#' @examples
cdf.laplace.regression = function(fit,data){
  coeff.hat=theta.hat[-1]
  sd.hat=theta.hat[1]
  y.hat = x%*%coeff.hat
  pnorm(y,mean=y.hat,sd=sd.hat)
}
