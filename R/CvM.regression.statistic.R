#' Cramer-von Mises statistic for regression
#'
#' @description
#'
#' @param y response variable
#' @param x explanatory variables
#' @param parameter estimates of regression parameters
#'
#' @return
#'
#' @seealso
#'
#' @name CvM.regression
#' @examples
#'
NULL

#' @export
#' @rdname CvM.regression
CvM.weibull.regression <- function(y,x,parameter=estimate.weibull.regression(y,x)){
  pp=length(parameter)
  p=pp-1
  beta = parameter[-pp]
  sigma = parameter[pp]
  shape=1/sigma
  scale = exp(x %*% beta)
  z <- pweibull(y,shape=shape,scale =scale)
  CvM(z)
}


# Helpers -----------------------------------------------------------------


cdf.normal.regression = function(y,x,theta.hat){
  #' Computes the probability integral transforms of the responses y
  #' in a normal linear regression model
  coeff.hat=theta.hat[-1]
  sd.hat=theta.hat[1]
  y.hat = x%*%coeff.hat
  pnorm(y,mean=y.hat,sd=sd.hat)
}


cdf.gamma.regression = function(y,x,thetahat,link="log"){
  n = length(y)
  pp = length(thetahat)
  coefs = thetahat[-pp]
  alpha=thetahat[pp]
  if( link == "log") invlink = exp
  if( link == "inverse" ) invlink = function(w) 1/w
  if( link == "identity" ) invlink = function(w) 1

  eta = x %*% as.matrix(coefs,nrow=pp-1)
  mu = invlink( eta )
  pgamma(y,shape=alpha,scale = mu/alpha)
}


cdf.laplace.regression = function(fit,data){
  coeff.hat=theta.hat[-1]
  sd.hat=theta.hat[1]
  y.hat = x%*%coeff.hat
  pnorm(y,mean=y.hat,sd=sd.hat)
}






