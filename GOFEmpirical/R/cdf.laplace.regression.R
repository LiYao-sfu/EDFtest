#' cdf.laplace.regression
#'
#' @param fit
#' @param data
#'
#' @return
#' @export
#'
#' @examples
cdf.laplace.regression = function(fit,data){
  #
  #  Computes the probability integral transforms of the responses y
  #   in a normal linear regression model
  #
  coeff.hat=theta.hat[-1]
  sd.hat=theta.hat[1]
  y.hat = x%*%coeff.hat
  pnorm(y,mean=yhat,sd=sd.hat)
}
