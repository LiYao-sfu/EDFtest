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
#' @export
#'
#' @examples
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
  Score=matrix(0,nrow=n,ncol=p+1)       # nrow = m? might be a typo
  Score[,1]=-1/sig.hat +(y-yhat)^2/sig.hat^3
  scaled.residual = sig.hat*(y-yhat)/yhat     # shape.hat might be sig.hat
  Score[,2:(p+1)]= (y-yhat)/sig.hat^2
  Score
}
