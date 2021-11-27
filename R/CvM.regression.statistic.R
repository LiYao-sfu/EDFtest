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
#' n = 500
#' p = 3
#' beta = c(1,2,3)
#' x = rnorm(n*(p))
#' #x = c(rep(1,n),x)
#' x = matrix(x,n,p)
#' # OLS regression
#' mean = x%*%(beta)
#' #apply the link to get the mean
#' #generate data with that mean,
#'
#' sd = 2
#' y =rnorm(n,mean=mean,sd=sd)
#' estimate.normal.regression(x=x,y=y,fit.intercept=TRUE)
#' CvM.normal.regression(x=x,y=y,fit.intercept=TRUE) # mean(cvm)=.06
#'
#'
NULL


#' @export
#' @rdname CvM.regression
CvM.normal.regression <- function(y,x,fit.intercept = TRUE){
  if(fit.intercept){fit = lm(y~x)}else{fit=lm(y~x-1)}
  xx = model.matrix(fit)
  beta = coef(fit)
  sigma = sqrt(mean(fit$residuals^2))
  mu = xx %*% beta
  z <- pnorm(y,mean=mu,sd=sigma)
  list(w = CvM(z),x.design=xx,betahat=beta,sigmahat=sigma)
}


#' @export
#' @rdname CvM.regression

CvM.gamma.regression <- function(y,x,fit,fit.intercept = TRUE,
                                 link="log"){
  if(missing(fit)){
    if (missing(x) || missing(y))stop("No fit is provided and one of x and y is missing")
    if(fit.intercept) {fit = glm(y~x, family = Gamma(link = link))}
    else{fit = glm(y~x-1, family = Gamma(link = link))}
  }
  xx = model.matrix(fit)
  beta = coef(fit)
  shape = 1/summary(fit)$dispersion

  if( link == "log" ) invlink = exp
  if( link == "inverse" ) invlink = function(w) 1/w
  if( link == "identity" ) invlink = function(w) w
  mu = invlink(xx %*% beta)
  scale = mu/shape
  z <- pgamma(y,shape=shape,scale = scale)
  list(w = CvM(z), x.design = xx,betahat=beta, shapehat=shape)
}


#' @export
#' @rdname CvM.regression
CvM.logistic.regression <- function(y,x,fit.intercept = TRUE){
  w = estimate.logistic.regression(y,x,fit.intercept=fit.intercept)
  theta=w$thetahat
  pp=length(theta)
  xx=w$x.design
  p=pp-1
  beta = parameter[-pp]
  sigma = parameter[pp]
  mu = x %*% beta
  z <- plogis(y,location=mu,scale=sigma)
  list(w=CvM(z),x.design=xx,betahat=beta,sigmahat=sigma)
}


#' @export
#' @rdname CvM.regression
CvM.laplace.regression <- function(y,x,fit.intercept = TRUE,
                                   parameter=estimate.laplace.regression(y,x)){
  pp=length(parameter)
  if(fit.intercept){
    x=cbind(1,x) # cbind(rep(1,dim(x)[1]),x)
  }
  p=pp-1
  beta = parameter[-pp]
  sigma = parameter[pp]
  mu = x %*% as.matrix(beta)
  z <- rmutil::plaplace(y,m=mu,s=sigma)
  CvM(z)
}


#' @export
#' @rdname CvM.regression
CvM.weibull.regression <- function(y,x,fit.intercept = TRUE,
                                   parameter=estimate.weibull.regression(y,x)){
  pp=length(parameter)
  if(fit.intercept){
    x=cbind(1,x) # cbind(rep(1,dim(x)[1]),x)
  }
  p=pp-1
  beta = parameter[-pp]
  sigma = parameter[pp]
  shape=1/sigma
  scale = exp(x %*% as.matrix(beta))
  z <- pweibull(y,shape=shape,scale =scale)
  CvM(z)
}


#' @export
#' @rdname CvM.regression
CvM.extremevalue.regression <- function(y,x,fit.intercept = TRUE,
                                        parameter=estimate.extremevalue.regression(y,x)){
  pp=length(parameter)
  if(fit.intercept){
    x=cbind(1,x) # cbind(rep(1,dim(x)[1]),x)
  }
  p=pp-1
  beta = parameter[-pp]
  sigma = parameter[pp]
  yy=(y-x %*% beta)/sigma
  z = exp(-exp(yy))
  CvM(z)
}


#' @export
#' @rdname CvM.regression
CvM.exp.regression <- function(y,x,fit,fit.intercept = TRUE,
                                 link="log"){
  if(missing(fit)){
    if (missing(x) || missing(y))stop("No fit is provided and one of x and y is missing")
    if(fit.intercept) {fit = glm(y~x, family = Gamma(link = link))}
    else{fit = glm(y~x-1, family = Gamma(link = link))}
  }
  xx = model.matrix(fit)
  beta = coef(fit)
  estimates = estimate.exp.regression(y,x,fit,
                                    fit.intercept = fit.intercept,
                                    link = link)
  beta=estimates$thetahat

  if( link == "log" ) invlink = exp
  if( link == "inverse" ) invlink = function(w) 1/w
  if( link == "identity" ) invlink = function(w) w
  mu = invlink(xx %*% beta)
  z <- pexp(y,rate = 1/mu)
  list(w = CvM(z), x.design = xx,betahat=beta)


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






