#' AD.logistic.eigen
#'
#' @param n number of eigenvalues
#'
#' @return
#' @export
AD.logistic.eigen = function(n){
  mean.asq.logistic = 1-(2*pi^2-3)/(2*(pi^2+3))  # from Maple
  M=AD.logistic.covmat(n)
  e=eigen(M)$values/n
  e*mean.asq.logistic/sum(e)
}
