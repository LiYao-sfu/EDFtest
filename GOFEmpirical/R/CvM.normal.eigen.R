#' CvM.normal.eigen
#'
#' @param n number of eigenvalues
#'
#' @return
CvM.normal.eigen = function(n){
  mean.wsq.normal=1/6 -7*sqrt(3)/(36*pi)
  M=CvM.normal.covmat(n)
  e=eigen(M)$values/n
  e*mean.wsq.normal/sum(e)
}
