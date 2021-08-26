#' AD.normal.eigen
#'
#' @param n number of eigenvalues
#'
#' @return
#' @export
AD.normal.eigen = function(n){
  M=AD.normal.covmat(n)
  e=eigen(M)$values/n
  e
}
