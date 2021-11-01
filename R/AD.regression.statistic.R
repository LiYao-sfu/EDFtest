#' Anderson-Darling statistic for regression
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
#' @name AD.regression
#' @examples
#'
NULL

#' @export
#' @rdname AD.regression
AD.extremevalue.regression <- function(y,x,parameter=estimate.extremevalue.regression(y,x)){
  pp=length(parameter)
  p=pp-1
  beta = parameter[-pp]
  sigma = parameter[pp]
  yy=(y-x %*% beta)/sigma
  z = exp(-exp(yy))
  AD(z)
}
