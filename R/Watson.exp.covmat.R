#' Watson.exp.covmat
#'
#' @param n number of eigenvalues
#'
#' @return
Watson.exp.covmat=function(n){
  s=1:n
  s=s/(n+1)
  (diag(n)-matrix(1/n,n,n)) %*% CvM.exp.covmat(n) %*% (diag(n)-matrix(1/n,n,n))
}
