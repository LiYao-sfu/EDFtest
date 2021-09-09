#' CvM.exp.eigen
#'
#' @param n number of eigenvalues
#'
#' @return
CvM.exp.eigen = function(n){
  mean = 5/54 # from Maple
  M=CvM.exp.covmat(n)
  e=eigen(M)$values/n
  e  * mean / sum(e)
}
