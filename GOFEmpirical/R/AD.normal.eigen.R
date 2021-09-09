#' AD.normal.eigen
#'
#' @param n number of eigenvalues
#'
#' @return
AD.normal.eigen = function(n){
  M=AD.normal.covmat(n)
  e=eigen(M)$values/n
  e
}
