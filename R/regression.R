#' Probability Integral Transforms for Normal Regression
#'
#' Computes the probability integral transforms of the responses y
#' in a normal linear regression model
#'
#' @param y response variable
#' @param x explanatory variables
#' @param theta.hat estimates of coefficients
#'
#' @return cdf.normal.regression gives probability integral transforms for y in a ordinary linear regression model
#'
#' @examples
cdf.normal.regression = function(y,x,theta.hat){
  coeff.hat=theta.hat[-1]
  sd.hat=theta.hat[1]
  y.hat = x%*%coeff.hat
  pnorm(y,mean=y.hat,sd=sd.hat)
}

#' Probability Integral Transforms for Gamma Regression
#'
#' Computes the probability integral transforms of the responses y
#' in a gamma log-linear regression model: log(E(y)) = x beta
#'
#' @param y response variable
#' @param x explanatory variables
#' @param theta.hat estimates of coefficients
#'
#' @return cdf.gamma.regression gives probability integral transforms for y in a gamma log-linear regression model
#'
#' @examples
cdf.gamma.regression = function(y,x,theta.hat){
  coeff.hat=theta.hat[-1]
  shape.hat=theta.hat[1]
  linearpredictor = x%*%coeff.hat
  yhat=exp(linearpredictor)
  pgamma(y,shape=shape.hat,scale=yhat/shape.hat)
}

#' Probability Integral Transforms for Laplace Regression
#'
#' Computes the probability integral transforms of the responses y
#' in a normal linear regression model
#'
#' @param y response variable
#' @param x explanatory variables
#' @param theta.hat estimates of coefficients
#'
#' @return cdf.laplace.regression gives probability integral transforms for y in a
#'
#' @examples
cdf.laplace.regression = function(fit,data){
  coeff.hat=theta.hat[-1]
  sd.hat=theta.hat[1]
  y.hat = x%*%coeff.hat
  pnorm(y,mean=y.hat,sd=sd.hat)
}

#' Fisher Information for Normal Regression Model
#'
#' The value returned is for sigma=1.  In general the FI must be divided by sigma^2.
#'
#' @param x matrix of covariates, normally x will contain an intercept term, that is a column of 1s
#'
fisher.information.normal=function(x){
  n=dim(x)[1]
  p=dim(x)[2]
  FI=matrix(0,nrow=p+1,ncol=p+1)
  FI[1,1]=1/2
  FI[2:(p+1),2:(p+1)]  = t(x)%*%x / n
  FI
}

#' Score Function for Simple Linear Regression Model
#'
#' computes an n by p+1 matrix containing the components of the score function
#' in a linear regression of y on the n x p covariate matrix x.
#' The ith row of the output corresponds to the ith data point.
#'
#' @param y the univariate response
#' @param x the design matrix
#' @param theta.hat MLEs of parameters
#'
#' @return Column 1 of the output is the sigma component of the score. Columns 2 to p+1
#' are the components of the score corresponding to the regression coefficients.
#'
#' @examples
#' m = matrix(rnorm(20),10,2)
#' X = cbind(rep(1,10),m)
#' y = rnorm(10)
#' fit = lm(y~X-1)
#' sigma = sqrt(mean(fit$residual^2))
#' theta.hat =  c(sigma,fit$coefficients)
#' score.normal.regression(y,X,theta.hat=theta.hat)
#' round(apply(.Last.value,2,sum))
score.normal.regression = function(y,x,theta.hat){
  #
  # As a check the column sums of the output should be 0 up to round off
  #  when the inputs labelled "theta.hat" are the MLEs of those parameters.
  #
  # But the sum of column 1 will be -p/sig.hat usually because sig.hat
  #  will be adjusted for degrees of freedom
  #
  # Should this be modified to expect a "fit" object from lm?
  #
  # If the model has an intercept this must be reflected in x having
  #  a column of 1s.
  #
  coeff.hat=theta.hat[-1]
  sig.hat=theta.hat[1]
  yhat = x%*%coeff.hat
  n=dim(x)[1]
  p=dim(x)[2]
  Score=matrix(0,nrow=n,ncol=p+1)
  Score[,1]=-1/sig.hat +(y-yhat)^2/sig.hat^3
  Score[,2:(p+1)]= x*(y-yhat)/sig.hat^2
  Score
}

#' Score Function for Gamma Regression Model
#'
#' Compute a p+1 by n matrix containing the components of the score function
#' in a Gamma regression of y on the n x p covariate matrix x.
#' The ith column of the output corresponds to the ith data point.
#'
#' @param y the univariate response
#' @param x the design matrix
#' @param theta.hat MLEs of parameters
#'
#' @return Row 1 of the output is the shape component of the score. Rows 2 to p+1
#' are the components of the score corresponding to the regression coefficients.
#'
#' @examples
#' set.seed(5)
#' n = 200  #(sample size)
#' p = 3
#' beta = c(3,4,5)
#' alpha = 3
#' x = rnorm(n*(p-1))
#' x = c(rep(1,n),x)
#' x = matrix(x,n,p)
#' scale=exp(x %*% beta) / alpha
#' y =rgamma(n,shape=alpha,scale=scale)
#' theta.hat = estimate.gamma.regression(x,y)
#'

score.gamma.regression = function(y,x,theta.hat){
  #
  # This computes a p+1 by n matrix containing the components of the
  #  score function in a Gamma regression of y on the n x p covariate matrix x
  # The ith column of the output corresponds to the ith data point.
  # Row 1 of the output is the shape component of the score.
  # Rows 2 to p+1 are the components of the score corresponding to the
  #  regression coefficients.
  #
  # As a check the row sums of the output should be 0 up to round off
  #  when the inputs labelled "?.hat" are the MLEs of those parameters.
  #
  coeff.hat=theta.hat[-1]
  shape.hat=theta.hat[1]
  linearpredictor = x%*%coeff.hat
  yhat=exp(linearpredictor)
  n=dim(x)[1]
  p=dim(x)[2]
  Score=matrix(0,nrow=p+1,ncol=n)
  Score[1,]=1+log(shape.hat)-digamma(shape.hat)-y/yhat +log(y/yhat)
  scaled.residual = shape.hat*(y-yhat)/yhat
  Score[2:(p+1),]= t(c(scaled.residual)*x)
  Score
}

