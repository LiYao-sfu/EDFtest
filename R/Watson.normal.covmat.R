#' Watson.normal.covmat
#'
#' @param n
#'
#' @return
#'
#' @examples
#'
Watson.normal.covmat=function(n){
  s=1:n
  s=s/(n+1)
  t=s
  CvM.normal.covmat(n) %*% (diag(n)-matrix(1/n,n,n))
  #(diag(n)-matrix(1/n,n,n)) %*% CvM.normal.covmat(n) %*% (diag(n)-matrix(1/n,n,n)) try this if not symmetric
  #row sum and col sum are zeros
}
