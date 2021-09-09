#' CvM.weibull.eigen
#'
#' @param n number of eigenvalues
#'
#' @return
CvM.weibull.eigen = function(n){
  mean.wsq.weibull=(1/54)-(4/9)*((log(3))^2-log(3)-1)/pi^2  # from Maple
  M=CvM.weibull.covmat(n)
  e=eigen(M)$values/n
  e*mean.wsq.weibull/sum(e)
}
