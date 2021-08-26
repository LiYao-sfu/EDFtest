#' AD.normal.covmat
#'
#' @param n number of eigenvalues
#'
#' @return
#' @export
AD.normal.covmat=function(n){
  s=1:n
  s=s/(n+1)
  t=s
  CvM.normal.covmat(n)/sqrt(outer(s*(1-s),t*(1-t)))
}
