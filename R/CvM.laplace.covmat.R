#' CvM.laplace.covmat
#'
#' @param n number of eigenvalues
#'
#' @return
CvM.laplace.covmat=function(n){
  Fisher.laplace = matrix(c(1,0,0,1),nrow=2)
  s=1:n
  s=s/(n+1)
  M1=outer(s,s,pmin)-outer(s,s)
  G1 = -s
  G1[s>0.5]=(s-1)[s>0.5]
  G2 = -s*log(2*s)
  G2[s>0.5] = ((1-s)*log(2*(1-s)))[s>0.5]
  M2 = cbind(G1,G2)
  M1-M2%*%solve(Fisher.laplace,t(M2))
}
