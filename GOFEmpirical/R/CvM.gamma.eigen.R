#' CvM.gamma.eigen
#'
#' @param n
#' @param shape
#'
#' @return
#' @export
#'
#' @examples
CvM.gamma.eigen = function(n,shape){
  # I don't have code for this yet
  # mean = 1/6 -?
  M=CvM.gamma.covmat(n)
  e=eigen(M)$values/n
  e # *mean/sum(e)
}
