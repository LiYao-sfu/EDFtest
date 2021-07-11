#' estimate.laplace.regression
#'
#' @param x
#' @param y
#' @param fit.intercept
#'
#' @return
#' @export
#'
#' @examples
estimate.laplace.regression=function(x,y,fit.intercept=TRUE){
  #
  # Just fits a standard linear model by least absolute deviations
  #  y is the univariate response and
  #  x is the design matrix which is
  #   NOT expected to have a column of 1s
  #  An intercept term is added by lm
  #
  # By default fits an intercept
  #
  data=data.frame(y=y,x=x)
  if(fit.intercept)fit=lad(y~x,data=data) else fit = lad(y~x-1,data=data)
  n = length(y)
  r = residuals(fit)
  coeff.hat = coefficients(fit)
  sigma.hat = sqrt(mean(abs(r)))
  c(sigma.hat,coeff.hat)
}
