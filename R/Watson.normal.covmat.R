#' Watson.normal.covmat
#'
#' @param n number of eigenvalues
#'
#' @return
Watson.normal.covmat=function(n){
  (diag(n)-matrix(1/n,n,n)) %*% CvM.normal.covmat(n) %*% (diag(n)-matrix(1/n,n,n))
}
