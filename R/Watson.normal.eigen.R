#' Watson.normal.eigen
#'
#' @param n number of eigenvalues
#'
#' @return
Watson.normal.eigen = function(n){
  M=Watson.normal.covmat(n)
  eigen(M)$values/n
}
