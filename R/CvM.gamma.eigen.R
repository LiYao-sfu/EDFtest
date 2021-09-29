#' CvM.gamma.eigen
#'
#' @param n number of eigenvalues
#' @param shape shape parameter of Gamma distribution
#'
#' @return
CvM.gamma.eigen = function(n,shape){
  # I don't have code for this yet
  # mean = 1/6 -?
  M=CvM.gamma.covmat(n,shape)
  e=eigen(M)$values/n
  e # *mean/sum(e)
}
