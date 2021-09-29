#' CvM.laplace.eigen
#'
#' @param n number of eigenvalues
#'
#' @return
CvM.laplace.eigen  = function(n){
  mean = 1/6 -1/12-1/54
  M=CvM.laplace.covmat(n)
  e=eigen(M)$values/n
  e * mean / sum(e)
}
