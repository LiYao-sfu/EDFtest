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
