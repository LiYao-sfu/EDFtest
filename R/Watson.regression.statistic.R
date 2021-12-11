#' Watson statistic for regression
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
#' @name Watson.regression
#' @examples
#'
NULL

#' @export
#' @rdname Watson.regression
Watson.normal.regression <- function(y,x,fit.intercept = TRUE){
  if(fit.intercept){fit = lm(y~x)}else{fit=lm(y~x-1)}
  xx = model.matrix(fit)
  beta = coef(fit)
  sigma = sqrt(mean(fit$residuals^2))
  mu = xx %*% beta
  z <- pnorm(y,mean=mu,sd=sigma)
  list(u=Watson(z),x.design=xx,betahat=beta,sigmahat=sigma)
}


#' @export
#' @rdname Watson.regression

Watson.gamma.regression <- function(y,x,fit,fit.intercept = TRUE,
                                 link="log"){
  w=estimate.gamma.regression(y,x,fit,fit.intercept = fit.intercept)
  xx = w$model.matrix
  theta =w$thetahat
  pp=length(theta)
  beta = theta[-pp]
  shape = theta[pp]

  if( link == "log" ) invlink = exp
  if( link == "inverse" ) invlink = function(w) 1/w
  if( link == "identity" ) invlink = function(w) w
  mu = invlink(xx %*% beta)
  scale = mu/shape
  z <- pgamma(y,shape=shape,scale = scale)
  list(u=Watson(z),x.design=xx,betahat=beta,shapehat=shape)
}


#' @export
#' @rdname Watson.regression
Watson.logistic.regression <- function(y,x,fit.intercept = TRUE){
  w = estimate.logistic.regression(y,x,fit.intercept=fit.intercept)
  theta=w$thetahat
  pp=length(theta)
  xx=w$model.matrix
  p=pp-1
  beta = theta[-pp]
  sigma = theta[pp]
  mu = xx %*% beta
  z <- plogis(y,location=mu,scale=sigma)
  list(u=Watson(z),x.design=xx,betahat=beta,sigmahat=sigma)
}


#' @export
#' @rdname Watson.regression
Watson.laplace.regression <- function(y,x,fit.intercept = TRUE){
  w = estimate.laplace.regression(y,x,fit.intercept=fit.intercept)
  theta=w$thetahat
  pp=length(theta)
  xx=w$model.matrix
  beta = theta[-pp]
  sigma = theta[pp]
  mu = xx %*% as.matrix(beta)
  #z <- rmutil::plaplace(y,m=mu,s=sigma)
  z = L1pack::plaplace(y,location = mu,scale = sigma)
  list(u=Watson(z),x.design=xx,betahat=beta,sigmahat=sigma)
}


#' @export
#' @rdname Watson.regression
Watson.weibull.regression <- function(y,x,fit.intercept = TRUE){
  w = estimate.weibull.regression(y,x,fit.intercept=fit.intercept)
  theta=w$thetahat
  pp=length(theta)
  xx=w$model.matrix
  beta = theta[-pp]
  shape = 1/theta[pp]
  scale = exp(xx %*% beta)
  z <- pweibull(y,shape=shape,scale =scale)
  list(u=Watson(z),x.design=xx,betahat=beta,shapehat=shape)
}


#' @export
#' @rdname Watson.regression
Watson.extremevalue.regression <- function(y,x,fit.intercept = TRUE){
  w = estimate.extremevalue.regression(y,x,fit.intercept=fit.intercept)
  theta=w$thetahat
  pp=length(theta)
  xx=w$model.matrix
  beta = theta[-pp]
  sigma = theta[pp]
  yy=(y-xx %*% beta)/sigma
  z = exp(-exp(yy))
  list(u=Watson(z),x.design=xx,betahat=beta,sigmahat=sigma)
}


#' @export
#' @rdname Watson.regression
Watson.exp.regression <- function(y,x,fit,fit.intercept = TRUE,
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
  list(u=Watson(z),x.design=xx,betahat=beta)
}
