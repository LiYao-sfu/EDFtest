#' estimate.normal.regression
#'
#' @param x
#' @param y
#' @param fit.intercept
#'
#' @return
#' @export
#'
#' @examples
estimate.normal.regression=function(x,y,fit.intercept=TRUE){
  #
  # Just fits a standard linear model by OLS
  #  y is the univariate response and
  #  x is the design matrix which is
  #   NOT expected to have a column of 1s
  #  An intercept term is added by lm
  #
  # By default fits an intercept
  #

  if(fit.intercept)fit=lm(y~x) else fit = lm(y~x-1)
  n = length(y)
  r = residuals(fit)
  coeff.hat = coefficients(fit)
  sigmahat = sqrt(mean(r^2)*n/(n-fit$df.residual))
  c(sigma.hat,coeff.hat)
}
