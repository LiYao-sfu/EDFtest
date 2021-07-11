#' CvM.laplace.eigen
#'
#' @param n
#'
#' @return
#' @export
#'
#' @examples
CvM.laplace.eigen  = function(n){
  mean = 1/6 -1/12-1/54
  M=CvM.laplace.covmat(n)
  e=eigen(M)$values/n
  e * mean / sum(e)
}
