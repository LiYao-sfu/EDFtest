#' Watson.laplace.covmat
#'
#' @param n number of eigenvalues
#'
#' @return
Watson.laplace.covmat=function(n){
  s=1:n
  s=s/(n+1)
  s2=sqrt(s*(1-s))
  (diag(n)-matrix(1/n,n,n)) %*% CvM.laplace.covmat(n) %*% (diag(n)-matrix(1/n,n,n))
}
