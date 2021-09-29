#' AD.exp.eigen
#'
#' @param n number of eigenvalues
#'
#' @return
AD.exp.eigen=function(n){
  mean = 0.595886194 # 3 - 2*zeta(3) from Maple
  s=1:n
  s=s/(n+1)
  M=AD.exp.covmat(n)
  e=eigen(M)$values/n
  e * mean / sum(e)
}
