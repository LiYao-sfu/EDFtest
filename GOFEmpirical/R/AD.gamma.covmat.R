#' AD.gamma.covmat
#'
#' @param n number of eigenvalues
#' @param shape
#'
#' @return
#' @export
AD.gamma.covmat=function(n,shape){
  s=1:n
  s=s/(n+1)
  CvM.gamma.covmat(n,shape=shape)/sqrt(outer(s*(1-s),s*(1-s)))
}
