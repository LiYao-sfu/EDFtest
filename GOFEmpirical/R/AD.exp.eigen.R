#' AD.exp.eigen
#'
#' @param n
#' @param shape
#'
#' @return
#' @export
#'
#' @examples
AD.exp.eigen=function(n,shape){
  #
  mean = 0.595886194 # 3 - 2*zeta(3) from Maple
  s=1:n
  s=s/(n+1)
  M=AD.exp.covmat(n)
  e=eigen(M)$values/n
  e * mean / sum(e)
}
