#' AD.laplace.covmat
#'
#' @param n number of eigenvalues
#'
#' @return
AD.laplace.covmat=function(n){
  s=1:n
  s=s/(n+1)
  s2=sqrt(s*(1-s))
  CvM.laplace.covmat(n)/outer(s2,s2)
}
