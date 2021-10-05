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
