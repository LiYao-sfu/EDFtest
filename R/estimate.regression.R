#' Estimate Regression Model Parameters
#'
#' @description
#' \code{estimate.normal.regression} fits a standard linear model by OLS (ordinary least square);
#' \code{estimate.laplace.regression} fits a standard linear model by LAD (least absolute deviations);
#' \code{estimate.gamma.regression} estimate the parameters in a Gamma regression model
#' which fits log(E(y))= x beta
#'
#' @param x A design matrix which is not expected to have a column of 1s.
#' @param y A univariate response.
#' @param fit.intercept Logical; if TRUE, an intercept term is added by lm.
#'
#' @return Estimated sigma and coefficients of a linear regression.
#'
#' @export
#' @examples
#' # OLS regression ------------------------------------
#' n = 500
#' p = 3
#' beta = c(1,2,3)
#' x = rnorm(n*(p-1))
#' x = c(rep(1,n),x)
#' x = matrix(x,n,p)
#' mean = x%*%(beta)
#' sd = 2
#' y =rnorm(n,mean=mean,sd=sd)
#' estimate.normal.regression(x,y)
#'
#' # LAD regression ------------------------------------
#' y =L1pack::rlaplace(n=n, location=mean, scale=sd)
#' estimate.laplace.regression(x,y,FALSE)
#'
#' # Gamma regression ----------------------------------
#' alpha = 3
#' scale=exp(x %*% beta) / alpha
#' y =rgamma(n=n,shape=alpha,scale=scale)
#' estimate.gamma.regression(x,y)
estimate.normal.regression=function(x,y,fit.intercept=TRUE){
  if(fit.intercept)fit=lm(y~x[,-1]) else fit = lm(y~x-1)
  n = length(y)
  r = residuals(fit)
  coeff.hat = coefficients(fit)
  sigma.hat = sqrt(mean(r^2)*n/(fit$df.residual))
  c(sigma.hat,coeff.hat)
}

#' @export
#' @rdname estimate.normal.regression
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

#' @export
#' @rdname estimate.normal.regression
estimate.gamma.regression=function(x,y){
  xx=x[,-1]
  fit=glm(y~xx,family=Gamma(link="log"))
  coeff.hat=coefficients(fit)
  yhat=fitted(fit)
  scaled.response=y/yhat
  target = mean(log(scaled.response))
  aold=1/var((y-yhat)/yhat)
  old.score=log(aold)-digamma(aold)+target
  old.score.derivative= 1/aold - trigamma(aold)
  anew <- aold - old.score/old.score.derivative
  if( anew < 0) anew <- aold/2
  while ( abs(anew-aold) > 1e-8){
    aold <- anew
    old.score = log(aold)-digamma(aold)+target
    old.score.derivative = 1/aold-trigamma(aold)
    anew <- aold - old.score/old.score.derivative
    if( anew < 0) anew <- aold/2
  }
  shape.hat=anew
  c(shape.hat, coeff.hat)
}
