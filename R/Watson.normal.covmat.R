#' Watson.normal.covmat
#'
#' @param n number of eigenvalues
#'
#' @return
Watson.normal.covmat=function(n){
  s=1:n
  s=s/(n+1)
  t=s
  #x=CvM.normal.covmat(n) %*% (diag(n)-matrix(1/n,n,n))
  (diag(n)-matrix(1/n,n,n)) %*% CvM.normal.covmat(n) %*% (diag(n)-matrix(1/n,n,n))
}
