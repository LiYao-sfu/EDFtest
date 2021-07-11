#' fisher.information.normal
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
fisher.information.normal=function(x){
  #
  # returns the Fisher Information per point
  # for a normal regression model: the mean is predicted
  # linearly from a matrix of covariates x
  # Normally x will contain an intercept term, that is,
  # a column of 1s
  #
  # The value returned is for sigma=1.  In general the FI must
  #  be divided by sigma^2.
  #
  n=dim(x)[1]
  p=dim(x)[2]
  FI=matrix(0,nrow=p+1,ncol=p+1)
  FI[1,1]=1/2
  FI[2:(p+1),2:(p+1)]  = t(x)%*%x / n
  FI
}
