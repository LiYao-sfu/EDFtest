#' AD.weibull.covmat
#'
#' @param n number of eigenvalues
#'
#' @return
#' @export
AD.weibull.covmat=function(n){
  s=1:n
  s=s/(n+1)
  t=s
  CvM.weibull.covmat(n)/sqrt(outer(s*(1-s),t*(1-t)))
}
