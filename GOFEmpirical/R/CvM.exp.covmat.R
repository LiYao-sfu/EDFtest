#' CvM.exp.covmat
#'
#' @param n
#'
#' @return
#' @export
#'
#' @examples
CvM.exp.covmat=function(n){
  FI = 1 # fisher.information.exp=1
  #
  s=1:n
  s=s/(n+1)
  M1=outer(s,s,pmin)-outer(s,s)
  G1 = -log(1-s)*(1-s)
  M2 = cbind(G1)
  M1-M2%*%solve(FI,t(M2))
}
