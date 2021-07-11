#' estimate.gamma.regression
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
estimate.gamma.regression=function(x,y){
  #
  # Estimates the parameters in a Gamma regression model
  #  which fits log(E(y))= x beta
  #
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
