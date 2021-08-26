#' AD.laplace.eigen
#'
#' @param n number of eigenvalues
#'
#' @return
#' @export
AD.laplace.eigen = function(n){
  mean = 1-.5351471355520514227 # from Maple
  M=AD.laplace.covmat(n)
  e=eigen(M)$values/n
  e  * mean / sum(e)
}
