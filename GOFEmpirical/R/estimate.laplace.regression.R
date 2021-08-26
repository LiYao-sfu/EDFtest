#' linear model by LAD
#'
#' Just fit a standard linear model by least absolute deviations
#'
#' @param x the design matrix which is not expected to have a column of 1s
#' @param y the univariate response
#' @param fit.intercept logical; if TRUE, an intercept term is added by lm
#'
#' @return estimated sigma and coefficients of linear regression by LAD.
#' @export
#'
#' @examples
#' library(LaplacesDemon)
#' set.seed(5)
#' n = 500  #(sample size)
#' p = 3
#' beta = c(1,2,3)
#' x = rnorm(n*(p-1))
#' x = c(rep(1,n),x)
#' x = matrix(x,n,p)
#' mean = x%*%(beta)
#' sd = 2
#' y =rlaplace(500, location=mean, scale=sd)
#' estimate.laplace.regression(x,y,FALSE)
estimate.laplace.regression=function(x,y,fit.intercept=TRUE){
  data=data.frame(y=y,x=x)
  require(L1pack)
  if(fit.intercept)fit=lad(y~x,data=data) else fit = lad(y~x-1,data=data)
  n = length(y)
  r = residuals(fit)
  coeff.hat = coefficients(fit)
  sigma.hat = fit$scale
  c(sigma.hat,coeff.hat)
}
