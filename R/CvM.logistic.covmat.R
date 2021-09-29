#' CvM.logistic.covmat
#'
#' @param n number of eigenvalues
#'
#' @return
#' @export
CvM.logistic.covmat=function(n){
  Fisher.logistic = matrix(c(1/3,0,0,(pi^2+3)/9),nrow=2)
  s=1:n
  s=s/(n+1)
  t=s
  M1=outer(s,t,pmin)-outer(s,t)
  G1 = s*(1-s)
  G2 = G1*(log(s/(1-s)))
  M2 = cbind(G1,G2)
  M1-M2%*%solve(Fisher.logistic,t(M2))
}
