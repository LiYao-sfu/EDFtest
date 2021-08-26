#' AD.exp.covmat
#'
#' @param n number of eigenvalues
#'
#' @return
#' @export
AD.exp.covmat=function(n){
  s=1:n
  s=s/(n+1)
  CvM.exp.covmat(n)/sqrt(outer(s*(1-s),s*(1-s)))
}
