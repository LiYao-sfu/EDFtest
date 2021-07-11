#' CvM.weibull.covmat
#'
#' @param n
#'
#' @return
#' @export
#'
#' @examples
CvM.weibull.covmat=function(n){
  Fisher.weibull = matrix(c(pi^2/6+(1+digamma(1))^2,-(1+digamma(1)),-(1+digamma(1)),1),nrow=2)
  s=1:n-0.5
  s=s/(n)
  t=s
  M1=outer(s,t,pmin)-outer(s,t)
  G1 = -log(1-s)*log(-log(1-s))*(1-s)
  G2 = log(1-s)*(1-s)
  M2 = cbind(G1,G2)
  M1-M2%*%solve(Fisher.weibull,t(M2))
  #sum(eigen(M1-M2%*%solve(Fisher.weibull,t(M2)))$values)/n
}
