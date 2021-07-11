#' CvM.logistic.eigen
#'
#' @param n
#'
#' @return
#' @export
#'
#' @examples
CvM.logistic.eigen = function(n){
  mean.wsq.logistic= 1/6 -(4*pi^2-9)/(20*(pi^2+3))  # from Maple
  M=CvM.logistic.covmat(n)
  e=eigen(M)$values/n
  e*mean.wsq.logistic/sum(e)
}
