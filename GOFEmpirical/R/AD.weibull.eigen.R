#' AD.weibull.eigen
#'
#' @param n number of eigenvalues
#'
#' @return
#' @export
AD.weibull.eigen = function(n){
  mean.asq.weibull = 0.3868394505  # from Maple
  M=AD.weibull.covmat(n)
  e=eigen(M)$values/n
  e*mean.asq.weibull/sum(e)
}
