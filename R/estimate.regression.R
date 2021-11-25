#' Estimate Regression Model Parameters
#'
#' @description
#' \code{estimate.normal.regression} fits a standard linear model by OLS (ordinary least square);
#' \code{estimate.laplace.regression} fits a standard linear model by LAD (least absolute deviations);
#' \code{estimate.gamma.regression} estimate the parameters in a Gamma regression model
#' which fits log(E(y))= x beta
#'
#' @param y A univariate response.
#' @param x A design matrix which is not expected to have a column of 1s.
#' @param fit.intercept Logical; if TRUE, an intercept term is added by lm.
#'
#' @return Estimated sigma and coefficients of a linear regression.
#'
#' @seealso
#'
#' @name estimate.regression
#' @examples
#' n = 500
#' p = 3
#' beta = c(1,2,3)
#' x = rnorm(n*(p-1))
#' x = c(rep(1,n),x)
#' x = matrix(x,n,p)
#'
#' # OLS regression
#' mean = x%*%(beta)
#' #apply the link to get the mean
#' #generate data with that mean,
#'
#' sd = 2
#' y =rnorm(n,mean=mean,sd=sd)
#' estimate.normal.regression(y=y,x=x)
#'
#' # LAD regression
#' library(L1pack)
#' y =L1pack::rlaplace(n=n, location=mean, scale=sd)
#' estimate.laplace.regression(y=y,x=x,FALSE)
#'
#' # Gamma regression
#' alpha = 3
#' scale=exp(x %*% beta) / alpha
#' y =rgamma(n=n,shape=alpha,scale=scale)
#' estimate.gamma.regression(y=y,x=x,intercept=FALSE)
#'
#' # Exponential regression
#' y =rexp(n=n,rate=1/exp(mean))
#' estimate.exp.regression(y=y,x=x)
#'
#' # Weibull regression
#' y = rweibull(n=n,shape=1,scale=exp(mean))
#' estimate.weibull.regression(y=y,x=x)
#'
#' # Extreme Value regression
#' y = log(rweibull(n=n,shape=1,scale=exp(mean)))
#' estimate.extremevalue.regression(y=y,x=x)
#'
#' # Logistic regression
#' y = rlogis(n=n,location=mean,scale=2)
#' estimate.logistic.regression(y=y,x=x)
#'
#'
#'
NULL


#' @export
#' @rdname estimate.regression
estimate.normal.regression=function(y,x,fit,fit.intercept=TRUE){
  if(fit.intercept)fit=lm(y~x) else fit = lm(y~x-1)
  n = length(y)
  r = residuals(fit)
  coeff.hat = coefficients(fit)
  sigma.hat = sqrt(mean(r^2)*n/(fit$df.residual))
  c(coeff.hat,sigma.hat)
}


#' @export
#' @rdname estimate.regression
estimate.gamma.regression = function(y,x,fit,fit.intercept=TRUE,link = "log"){
  #
  #  This function uses glm to get initial values for maximum
  #   likelihood fits of a model in which link(E(y)) =x %*% coefs
  #  We will test that the $y$ have gamma distributions with mean
  #    invlink(x %*% coefs)
  #  and shape alpha.  GLM gives us initial estimates of coefs and of
  #  the dispersion which is 1/alpha.
  #  Then optim is used to find the mle of coefs and of alpha
  #
  #  We allow the three link functions used by glm for the Gamma family
  #  So far we don't check that one of the possibilities is actually called.
  #  The code will crash if not.
  #

  if(missing(fit)){
    if (missing(x) || missing(y))stop("No fit is provided and one of x and y is missing")
    if(fit.intercept) {fit = glm(y~x, family = Gamma(link = link))}
    else{fit = glm(y~x-1, family = Gamma(link = link))}
  }
  xx = model.matrix(fit)
  betastart = coef(fit)
  shapestart = 1/summary(fit)$dispersion
  thetastart = c(betastart,shapestart)
  cat("Initial values ",thetastart,"\n")
  if( link == "log" ){
    invlink = exp
    weight = function(w) 1
    logh2der = function(w) 0
  }
  if( link == "inverse" ){
    invlink = function(w) 1/w
    weight = function(w) -1/w
    logh2der = function(w) 1/w^2
  }
  if( link == "identity" ){
    invlink = function(w) w
    weight = 1/w
    logh2der = function(w) -1/w^2
  }
  ellneg = function(th,xx,y){
    n = length(y)
    pp = length(th)
    coefs = th[-pp]
    alpha = th[pp]
    mu = invlink(xx %*% coefs)
    ell = alpha*log(y)+alpha*log(alpha/mu) -alpha*y/mu -lgamma(alpha)
    -sum(ell)
  }
  score = function(th,xx,y){
    n = length(y)
    pp = length(th)
    p = pp - 1
    coefs = th[-pp]
    alpha = th[pp]
    eta = xx %*% coefs
    mu = invlink(eta)
    wt = weight(eta)
    sc.alpha = log(y/mu) -y/mu -digamma(alpha)+log(alpha) +1
    sc.mu = xx*rep(wt*alpha*(y/mu -1),p)
    cbind(sc.mu,sc.alpha)
  }
  hessian = function(th,xx,y){
    n = length(y)
    pp = length(th)
    p = pp - 1
    coefs = th[-pp]
    alpha=th[pp]
    eta = xx %*% coefs
    mu = invlink(eta)
    wt = weight(eta)
    wt2 = logh2der(eta)
    sc.alpha = log(y/mu) -y/mu -digamma(alpha)+log(alpha) +1
    H = array(0,dim=c(n,pp,pp))
    h1.mu.mu = (y/mu-1)*wt2
    h2.mu.mu = -y/mu*wt^2
    h.mu.mu = alpha*(h1.mu.mu+h2.mu.mu)
    h.al.al = 1/alpha-trigamma(alpha)
    h.al.mu = (y/mu-1)*wt
    h.m.m = array(0,dim=c(n,p,p))
    for( i in 1:n){
      H[i,1:p,1:p] = outer(h.mu.mu[i]*xx[i,],xx[i,])
      H[i,1:p,pp] = h.al.mu[i]*xx[i,]
      H[i,pp,1:p] = h.al.mu[i]*xx[i,]
    }
    H[,pp,pp] = h.al.al + sc.alpha
    H
  }
  f = function(theta,predictor,response){
    ellneg(theta,predictor,response)
  }
  D1 = function(theta,predictor,response){
    -apply(score(theta,predictor,response),2,sum)
  }
  D2 = function(theta,predictor,response){
    -apply(hessian(theta,predictor,response),c(2,3),sum)
  }
  Marq = marqLevAlg::marqLevAlg(b=thetastart,fn=f,gr=D1,hess=D2,
                                 #print.info=TRUE,
                                 minimize = TRUE,maxiter=100,
                                 response = y,predictor=xx)
  thetahat = Marq$b
 # w = optim(par = thetastart, ellneg, xx=xx,y=y,invlink=invlink)
 # thetahat = w$par
 # list(thetahat = thetahat, model.matrix = xx,optim.output = w)
  list(thetahat = thetahat, model.matrix = xx,Marq.output = Marq)
}


#' @export
#' @rdname estimate.regression
estimate.exp.regression = function(y,x,fit,fit.intercept=TRUE,link = "log"){
  #
  #  This function uses glm to get initial values for maximum
  #   likelihood fits of a model in which link(E(y)) =x %*% coefs
  #  We will test that the $y$ have exponential distributions with mean
  #    invlink(x %*% coefs).
  #    GLM gives us initial estimates of coefs
  #  The dispersion  is 1.
  #  Then Marquardt-Levenberg is used to ensure we have found find the mle.
  #
  #  We allow the three link functions used by glm for the Gamma family
  #  So far we don't check that one of the possibilities is actually called.
  #  The code will crash if not.
  #
  if(missing(fit)){
    if(missing(x) || missing(y)) stop("No fit is provided and one of x and y is missing")
    if(fit.intercept) {fit = glm(y~x, family = Gamma(link = link))}
    else{fit = glm(y~x-1, family = Gamma(link = link))}
  }
  xx = model.matrix(fit)
  betastart = coef(fit)
  if( link == "log" ){
    invlink = exp
    weight = function(w) 1
    logh2der = function(w) 0
  }
  if( link == "inverse" ){
    invlink = function(w) 1/w
    weight = function(w) -1/w
    logh2der = function(w) 1/w^2
  }
  if( link == "identity" ){
    invlink = function(w) w
    weight = 1/w
    logh2der = function(w) - 1/w^2
  }
  ellneg = function(th,xx,y){
    n = length(y)
    pp = length(th)
    coefs = th
    mu = invlink(xx %*% coefs)
    ell = -log(mu) -y/mu
    -sum(ell)
  }
  score = function(th,xx,y){
    n = length(y)
    p = length(th)
    coefs = th
    eta = xx %*% coefs
    mu = invlink(eta)
    wt = weight(eta)
    sc.mu = xx*rep(wt*(-1 +y/mu),p)
    sc.mu
  }

  hessian = function(th,xx,y){
    n = length(y)
    p = length(th)
    coefs = th
    eta = xx %*% coefs
    mu = invlink(eta)
    wt = weight(eta)
    wt2 = logh2der(eta)
    H = array(0,dim=c(n,p,p))
    h1.mu.mu = (y/mu-1)*wt2
    h2.mu.mu = -y/mu*wt^2
    h.mu.mu = h1.mu.mu+h2.mu.mu
    for( i in 1:n){
      H[i,1:p,1:p] = outer(h.mu.mu[i]*wt[i]^2*xx[i,],xx[i,])
    }
    H
  }

  f = function(theta,predictor,response){
    ellneg(theta,predictor,response)
  }
  D1 = function(theta,predictor,response){
    -apply(score(theta,predictor,response),2,sum)
  }
  D2 = function(theta,predictor,response){
    -apply(hessian(theta,predictor,response),c(2,3),sum)
  }

  Marq = marqLevAlg::marqLevAlg(b=betastart,fn=f,gr=D1,hess=D2,
                                #print.info=TRUE,
                                minimize = TRUE,maxiter=100,
                                response = y,predictor=xx)
  thetahat = Marq$b
  list(thetahat = thetahat, model.matrix = xx,Marq.output = Marq)
}


#' @export
#' @rdname estimate.regression
estimate.logistic.regression <- function(x,y,fit,fit.intercept=TRUE,detail=FALSE){
  #
  # Use the Marquardt-Levenberg algorithm to fit a logistic regression
  #  model in which the log of the response is predicted by x
  # The scale parameter is taken to be x^\top beta .
  # To get initial values we do linear regression remembering
  #
  #
  if(fit.intercept)fit=lm(y~x) else fit = lm(y~x-1)
  beta.start = coef(fit)
  sigma.start =  summary(fit)$sigma*sqrt(3/pi^2)

  theta = c(beta.start,sigma.start)
  #  print(theta)
  f = function(theta,response,predictor){
    -ell.logistic.regression(response,predictor,theta)
  }
  D1 = function(theta,response,predictor){
    -apply(score.logistic.regression(response,predictor,theta),2,sum)
  }
  D2 = function(theta,response,predictor){
    -apply(hessianarray.logistic.regression(response,predictor,theta),c(2,3),sum)
  }
  # cat("Initial Log Likelihood ", f(theta,log(y),x),"\n")
  # cat("Initial Score ", D1(theta,log(y),x),"\n")
  # cat("Initial Hessian\n")
  # print(D2(theta,log(y),x))
  # print(eigen(D2(theta,log(y),x))$values)
  Marq = marqLevAlg::marqLevAlg(b=theta,fn=f,gr=D1,hess=D2,
                                #print.info=TRUE,
                                minimize = TRUE,maxiter=100,
                                response = y,predictor=x)
  thetahat = Marq$b
  # print(Marq)
  # print(rbind(thetahat,theta))
  if(detail) return(list(thetahat = thetahat, Marq = Marq))
  thetahat
}


#' @export
#' @rdname estimate.regression
estimate.laplace.regression=function(y,x,fit.intercept=TRUE){
  data=data.frame(y=y,x=x)
  if(fit.intercept)fit=lad(y~x,data=data) else fit = lad(y~x-1,data=data)
  n = length(y)
  r = residuals(fit)
  coeff.hat = coefficients(fit)
  sigma.hat = fit$scale
  c(coeff.hat,sigma.hat)
}


#' @export
#' @rdname estimate.regression
estimate.weibull.regression <- function(y,x,fit.intercept=TRUE,detail=FALSE){
  #
  # Use the Marquardt-Levenberg algorithm to fit a weibull regression
  #  model in which the log of the scale parameter for the response is predicted by x
  # That is, the scale parameter is taken to be log(x^\top beta).
  # To get initial values we do linear regression remembering
  #  that log(y) -x^\top beta has mean -gamma and variance pi^2/6 times the shape parameter
  #  where gamma is Euler's constant.
  #
  w=log(y)-digamma(1)
  if(fit.intercept){
    #
    # The design matrix x is assumed not to contain a column
    # of 1s and an intercept is added to the fit
    #
    fit = lm(w~x)
    beta.start = coef(fit)[-1]
    int.start = coef(fit)[1]+digamma(1)
    sigma.start =  summary(fit)$sigma*6/pi^2
    theta = c(int.start,beta.start,sigma.start)
    xx = cbind(rep(1,n),x)
    } else{
    #
    # It is assumed that a column of 1s is included in x or that
    # for some reason a column of 1s is not desired.
    #
    fit = lm(w~x-1)
    beta.start = coef(fit)
    sigma.start =  summary(fit)$sigma*6/pi^2
    theta = c(beta.start,sigma.start)
    xx=x
  }

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
                                minimize = TRUE,maxiter=500,response = log(y),predictor=xx)
  thetahat = Marq$b
  if(detail) return(list(thetahat = thetahat, Marq = Marq))
  thetahat
}


#' @export
#' @rdname estimate.regression
estimate.extremevalue.regression <- function(x,y,fit.intercept=TRUE,detail=FALSE){
  #
  # Use the Marquardt-Levenberg algorithm to fit an Extreme value regression
  #  model in which the response is predicted by x
  # The scale parameter is taken to be x^\top beta .
  # To get initial values we do linear regression remembering
  #  that y -x^\top beta has mean -gamma and variance pi^2/6
  #  where gamma is Euler's constant.
  #
  w=y-digamma(1)
  if(fit.intercept){
    fit = lm(w~x)
    beta.start = coef(fit)[-1]
    int.start = coef(fit)[1]+digamma(1)
    sigma.start =  summary(fit)$sigma*6/pi^2
    xx=cbind(rep(1,length(y)),x)
    theta = c(int.start,beta.start,sigma.start)
  }else{
    fit = lm(w~x-1)
    beta.start = coef(fit)
    sigma.start =  summary(fit)$sigma*6/pi^2
    xx=x
    theta = c(beta.start,sigma.start)
   }

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
                                minimize = TRUE,maxiter=500,response = y,predictor=xx)
  thetahat = Marq$b
  if(detail) return(list(thetahat = thetahat, Marq = Marq))
  thetahat
}


# Helpers -------------------------------------


score.logistic.regression = function(y,x,theta){
  pp=length(theta)
  n=length(y)
  p=pp-1
  beta = theta[1:p]
  sigma = theta[pp]
  mu = x %*% beta
  r = (y-mu) / sigma
  er = exp(r)
  ua = x*rep((2*er/(1+er)-1) / sigma, p)
  us = ( (r-1)*er -r -1) / (sigma*(1+er))
  cbind(ua,us)
}

hessianarray.logistic.regression = function(y,x,theta){
  pp=length(theta)
  n=length(y)
  p=pp-1
  beta = theta[1:p]
  sigma = theta[pp]
  mu = x %*% beta
  r = (y-mu) / sigma
  er = exp(r)
  umbit = -2*er/(1+er)^2
  H = array(0,dim=c(n,pp,pp))
  umm = array(0,dim=c(n,p,p))
  for( i in 1:n) umm[i,,]= outer(umbit[i]*x[i,]/sigma^2,x[i,])
  ums = x*rep((1-2*r*er-er^2)/(sigma^2*(1+er)^2),p)
  uss = (1+2*r-(2*r-1)*er^2-2*(r^2-1)*er)/(sigma^2*(1+er)^2)
  H[,1:p,1:p] = umm
  H[,pp,1:p] = ums
  H[,1:p,pp] = ums
  H[,pp,pp]  = uss
  H
}

ell.logistic.regression = function(y,x,theta){
  pp=length(theta)
  n=length(y)
  p=pp-1
  beta = theta[1:p]
  sigma = theta[pp]
  mu = x %*% beta

  r = (y-mu) / sigma
  er = exp(r)
  ell =-n*log(sigma) + sum(r-2*log(1+er))
  ell
}


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





