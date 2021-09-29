#' linear model by OLS
#'
#' Just fit a standard linear model by OLS (ordinary least square)
#'
#' @param x the design matrix which is not expected to have a column of 1s
#' @param y the univariate response
#' @param fit.intercept logical; if TRUE, an intercept term is added by lm
#'
#' @return estimated sigma and coefficients of linear regression by OLS.
#' @export
#'
#' @examples
#' set.seed(5)
#' n = 100
#' p = 3
#' beta = c(1,2,3)
#' x = rnorm(n*(p-1))
#' x = c(rep(1,n),x)
#' x = matrix(x,n,p)
#' mean = x%*%(beta)
#' sd = 2
#' y =rnorm(n,mean=mean,sd=sd)
#' estimate.normal.regression(x,y)
estimate.normal.regression=function(x,y,fit.intercept=TRUE){
  if(fit.intercept)fit=lm(y~x[,-1]) else fit = lm(y~x-1)
  n = length(y)
  r = residuals(fit)
  coeff.hat = coefficients(fit)
  sigma.hat = sqrt(mean(r^2)*n/(fit$df.residual))
  c(sigma.hat,coeff.hat)
}
