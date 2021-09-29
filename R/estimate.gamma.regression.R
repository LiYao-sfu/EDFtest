#' Gamma regression model
#'
#' Estimate the parameters in a Gamma regression model which fits log(E(y))= x beta
#'
#' @param x the design matrix
#' @param y the univariate response
#'
#' @return estimated sigma and coefficients of gamma regression model.
#' @export
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
#' estimate.gamma.regression(x,y)
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
