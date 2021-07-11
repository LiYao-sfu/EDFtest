#' gof.sandwich
#'
#' @param y
#' @param x
#' @param Fdist
#' @param thetahat
#' @param Score
#' @param m
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
gof.sandwich=function(y,x=NULL,Fdist,thetahat,Score,m=max(n,100),...){
  require(CompQuadForm)
  #
  #  This function tests the hypothesis that data y come from
  #  distribution Fdist with unknown parameter values theta
  #
  #  Estimates of theta must be provided in thetahat.
  #
  #  It uses a large sample approximation to the limit distribution
  #  based on the use of the score function components
  #  to estimate the Fisher information and the limiting covariance
  #  function of the empirical process.
  #
  #  The estimates thetahat should be roots of the likelihood equations.
  #
  #  Input:
  #         y         data -- should be a numerical vector
  #         Fdist     the hypothesized cdf -- must return a vector of
  #                     probability integral transforms of length = length(y)
  #         thetahat  parameter estimates -- the mles
  #         Score     returns components of the score function
  #                     an n by p matrix with entries partial log f(x_i,\theta)/ partial theta_j
  #         m         Eigenvalues are extracted for an m by m grid of the covariance ftn
  #         Other inputs passed to Fdist and Score when needed.
  #
  n = length(y)                # Sample size
  p=length(thetahat)           # Number of parameters
  if(is.null(x)){
    pit = Fdist(y,thetahat)   # Computes the probability integral transforms
  }else{
    pit=Fdist(y,x,thetahat)
  }
  if(is.null(x)){
    u =t(Score(y,thetahat))         # Components of the score: n by p matrix
  }else{
    u =t(Score(y,x,thetahat))         # Components of the score: n by p matrix
  }
  Fisher = t(u)%*% u / n       # Estimate of Fisher information in 1 point.
  s=(1:m)/(m+1)                # Grid on which to compute covariance matrix.
  #
  ind=function(x,y){
    one = rep(1,length(x))
    one[x>y] = 0
    one
  }
  #
  Dfb = outer(pit,s,ind)
  #
  #   Dfb is a matrix whose entry i,j is 1 if pit[i] > s[j]
  #       it is n by m
  #
  Df=t(Dfb) %*% u/n  # This is an m by p matrix.
  m1 = outer(s,s,pmin)
  m2 = outer(s,s,"*")
  Sigma.W =  m1-m2-Df%*%solve(Fisher,t(Df))
  Sigma.A =  Sigma.W/sqrt(outer(s*(1-s),s*(1-s),"*"))
  J = diag(m) - matrix(1/m,m,m)
  Sigma.U =  J %*% Sigma.W %*% J
  Evals.W = eigen(Sigma.W)$values/m
  Evals.A = eigen(Sigma.A)$values/m
  Evals.U = eigen(Sigma.U)$values/m
  stat=gof.statistics.only(pit)
  P.W = imhof(stat$W,Evals.W)$Qq
  P.A = imhof(stat$A,Evals.A)$Qq
  P.U = imhof(stat$U,Evals.U)$Qq
  #
  list(W2=list(W2=stat$W,P=P.W),A2=list(A2=stat$A,P=P.A),U2=list(U2=stat$U,P=P.U))
}
