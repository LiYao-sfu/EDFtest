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
#'
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


#'
#' @rdname estimate.normal.regression
estimate.gamma.regression = function(fit,x,y,link = "log"){
  #
  #  This function uses glm to get initial values for maximum
  #   likelihood fits of a model in which link(E(y)) =x %*% coefs
  #  We will test that the $y$ have gamma distributions with mean
  #    invlink(x %*% coefs)
  #  and shape alpha.  GLM gives us initial estimates of coefs and of
  #  the dispersion which is 1/alpha.
  #  Then optim is used to find the mle.
  #
  #  We allow the three link functions used by glm for the Gamma family
  #  So far we don't check that one of the possibilities is actually called.
  #  The code will crash if not.
  #
  if(missing(fit)){
    if (missing(x) || missing(y))stop("No fit is provided and one of x and y is missing")
    fit = glm(y~x, family = Gamma(link = link))
  }
  xx = model.matrix(fit)
  betastart = coef(fit)
  sigmastart = summary(fit)$dispersion
  thetastart = c(betastart,sigmastart)
  cat("Initial values ",thetastart,"\n")
  if( link == "log" ) invlink = exp
  if( link == "inverse" ) invlink = function(w) 1/w
  if( link == "identity" ) invlink = function(w) w
  ell = function(th,xx,y,invlink = invlink){
    n = length(y)
    pp = length(th)
    coefs = th[-pp]
    alpha=1/th[pp]
    mu = invlink(xx %*% coefs)
    ell = alpha*log(y)+alpha*log(alpha/mu) -alpha*y/mu -lgamma(alpha)
    -sum(ell)
  }
  w = optim(par = thetastart, ell, xx=xx,y=y,invlink=invlink)
  thetahat = w$par
  pp = length(thetahat)
  thetahat[pp] = 1/thetahat[pp]
  list(thetahat = thetahat, model.matrix = xx,optim.output = w)
}


#'
#' @rdname estimate.normal.regression
estimate.laplace.regression=function(x,y,fit.intercept=TRUE){
  data=data.frame(y=y,x=x)
  if(fit.intercept)fit=lad(y~x,data=data) else fit = lad(y~x-1,data=data)
  n = length(y)
  r = residuals(fit)
  coeff.hat = coefficients(fit)
  sigma.hat = fit$scale
  c(sigma.hat,coeff.hat)
}


estimate.weibull.regression <- function(y,x,detail=FALSE){
  #
  # Use the Marquardt-Levenberg algorithm to fit a weibull regression
  #  model in which the log of the response is predicted by x
  # The scale parameter is taken to be x^\top beta .
  # To get initial values we do linear regression remembering
  #  that log(y) -x^\top beta has mean -gamma and variance pi^2/6
  #  where gamma is Euler's constant.
  #
  w=log(y)-digamma(1)
  fit = lm(w~x-1)
  beta.start = coef(fit)[-1]
  int.start = coef(fit)[1]+digamma(1)
  sigma.start =  summary(fit)$sigma*6/pi^2

  theta = c(int.start,beta.start,sigma.start)

  f = function(theta,response,predictor){
    -ell.extremevalue.regression(response,predictor,theta)
  }
  D1 = function(theta,response,predictor){
    -apply(score.extremevalue.regression(response,predictor,theta),2,sum)
  }
  D2 = function(theta,response,predictor){
    -apply(hessianarray.extremevalue.regression(response,predictor,theta),c(2,3),sum)
  }
  Marq = marqLevAlg::marqLevAlg(b=theta,fn=f,gr=D1,hess=D2,
                                epsa=0.001,#print.info=TRUE,
                                minimize = TRUE,maxiter=500,response = log(y),predictor=x)
  thetahat = Marq$b
  if(detail) return(list(thetahat = thetahat, Marq = Marq))
  thetahat
}


estimate.extremevalue.regression <- function(y,x,detail=FALSE){
  #
  # Use the Marquardt-Levenberg algorithm to fit an Extreme value regression
  #  model in which the response is predicted by x
  # The scale parameter is taken to be x^\top beta .
  # To get initial values we do linear regression remembering
  #  that y -x^\top beta has mean -gamma and variance pi^2/6
  #  where gamma is Euler's constant.
  #
  w=y-digamma(1)
  fit = lm(w~x-1)
  beta.start = coef(fit)[-1]
  int.start = coef(fit)[1]+digamma(1)
  sigma.start =  summary(fit)$sigma*6/pi^2

  theta = c(int.start,beta.start,sigma.start)

  f = function(theta,response,predictor){
    -ell.extremevalue.regression(response,predictor,theta)
  }
  D1 = function(theta,response,predictor){
    -apply(score.extremevalue.regression(response,predictor,theta),2,sum)
  }
  D2 = function(theta,response,predictor){
    -apply(hessianarray.extremevalue.regression(response,predictor,theta),c(2,3),sum)
  }
  Marq = marqLevAlg::marqLevAlg(b=theta,fn=f,gr=D1,hess=D2,
                                epsa=0.001,#print.info=TRUE,
                                minimize = TRUE,maxiter=500,response = y,predictor=x)
  thetahat = Marq$b
  if(detail) return(list(thetahat = thetahat, Marq = Marq))
  thetahat
}

# Helpers -------------------------------------


score.extremevalue.regression = function(y,x,theta){
  pp=length(theta)
  n=length(y)
  p=pp-1
  beta = theta[1:p]
  sigma = theta[pp]
  mu = x %*% beta
  r = (y-mu) / sigma
  ua = x*rep((exp(r)-1) / sigma, p)
  us = ( r*exp(r) -r -1) / sigma
  cbind(ua,us)
}

hessianarray.extremevalue.regression = function(y,x,theta){
  pp=length(theta)
  n=length(y)
  p=pp-1
  beta = theta[1:p]
  sigma = theta[pp]
  mu = x %*% beta
  r = (y-mu) / sigma
  er = exp(r)
  umbit= -er/ sigma^2
  H = array(0,dim=c(n,pp,pp))
  umm = array(0,dim=c(n,p,p))
  for( i in 1:n) umm[i,,]= outer(umbit[i]*x[i,],x[i,])
  ums = -x*rep( (r*er +er -1) / sigma^2,p)
  uss = (-(r^2)*er -2*r*er+2*r + 1)/sigma^2
  H[,1:p,1:p] = umm
  H[,pp,1:p] = ums
  H[,1:p,pp] = ums
  H[,pp,pp]  = uss
  H
}

ell.extremevalue.regression = function(y,x,theta){
  pp=length(theta)
  n=length(y)
  p=pp-1
  beta = theta[1:p]
  sigma = theta[pp]
  mu = x %*% beta
  r = (y-mu) / sigma
  er = exp(r)
  # print(-n*log(sigma) + sum(r-er))
  (-n*log(sigma) + sum(r-er))
}





