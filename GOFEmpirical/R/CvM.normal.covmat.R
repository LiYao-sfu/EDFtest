#' CvM.normal.covmat
#'
#' @param n number of eigenvalues
#'
#' @return
CvM.normal.covmat=function(n){
  Fisher.normal = matrix(c(1,0,0,2),nrow=2)
  s=1:n
  s=s/(n+1)
  M1=outer(s,s,pmin)-outer(s,s)
  x = qnorm(s)
  G1 = dnorm(x)
  G2 = -x*G1
  M2 = cbind(G1,G2)
  M1-M2%*%solve(Fisher.normal,t(M2))
}
