#' CvM.gamma.covmat
#'
#' @param n number of eigenvalues
#' @param shape
#'
#' @return
#' @export
CvM.gamma.covmat=function(n,shape){
  fisher.information.gamma=function(shape.hat){
    #
    # returns the estimated Fisher Information per point
    # for a gamma regression model in which the log mean is predicted
    # linearly from a matrix of covariates x
    # Normally x will contain an intercept term
    #
    FI=matrix(0,nrow=2,ncol=2)
    FI[1,1]=trigamma(shape.hat)
    FI[1,2]= 1
    FI[2,1]= 1
    FI[2,2]=shape.hat
    FI
  }
  g = gamma(shape)
  dg = digamma(shape)
  FI = fisher.information.gamma(shape.hat = shape)
  s=1:n
  s=s/(n+1)
  M1=outer(s,s,pmin)-outer(s,s)
  G1 = s*0
  Q = qgamma(s,shape=shape)
  D = dgamma(Q,shape=shape)
  G2 = - Q * D
  g1.integrand = function(x,shape){log(x)*x^(shape-1)*exp(-x)}
  for(i in 1:n){
    G1[i] = integrate(g1.integrand,0,Q[i],shape=shape)$value/g -s[i]*dg
  }
  M2 = cbind(G1,G2)
  M1-M2%*%solve(FI,t(M2))
}
