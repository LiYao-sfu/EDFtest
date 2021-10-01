#' Watson.gamma.covmat
#'
#' @param n number of eigenvalues
#' @param shape shape parameter for Gamma distribution
#'
#' @return
Watson.gamma.covmat=function(n,shape){
  s=1:n
  s=s/(n+1)
  (diag(n)-matrix(1/n,n,n)) %*% CvM.gamma.covmat(n,shape=shape) %*% (diag(n)-matrix(1/n,n,n))
}
