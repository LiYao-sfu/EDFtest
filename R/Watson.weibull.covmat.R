#' Watson.weibull.covmat
#'
#' @param n number of eigenvalues
#'
#' @return
Watson.weibull.covmat=function(n){
  s=1:n
  s=s/(n+1)
  t=s
  (diag(n)-matrix(1/n,n,n)) %*% CvM.weibull.covmat(n) %*% (diag(n)-matrix(1/n,n,n))
}
