#' AD.gamma.eigen
#'
#' @param n number of eigenvalues
#' @param shape shape parameter for Gamma distribution
#'
#' @return
AD.gamma.eigen=function(n,shape){
  # I don't have code for this yet
  # mean = 1 -?
  s=1:n
  s=s/(n+1)
  M=AD.gamma.covmat(n,shape)
  e=eigen(M)$values/n
  e # *mean/sum(e)
}
