#' EDF Goodness-of-Fit tests for General Distributions using Sandwich Estimation of Covariance Function
#'
#' This function tests the hypothesis that data y come from
#' distribution Fdist with unknown parameter values theta
#'
#' Estimates of theta must be provided in thetahat.
#'
#' It uses a large sample approximation to the limit distribution
#' based on the use of the score function components
#' to estimate the Fisher information and the limiting covariance
#' function of the empirical process.
#'
#' The estimates thetahat should be roots of the likelihood equations.
#'
#' @param y data -- should be a numerical vector: sample or response of regression problem
#' @param x matrix of covariates
#' @param Fdist user supplied function to compute probability integral transform of y
#' @param thetahat parameter estimates by mle
#' @param Score user supplied function to compute3 components of the score function an n by p matrix with entries
#' partial log f(y_i,\theta)/ partial theta_j
#' @param m Eigenvalues are extracted for an m by m grid of the covariance ftn
#' @param ... other inputs passed to Fdist and Score when needed.
#'
#' @return
#' @export
#'
#' @examples
#'
gof.sandwich=function(y,x=NULL,Fdist,thetahat,Score,m=max(n,100),...){
  require(CompQuadForm)

  n = length(y)                   # Sample size
  p=length(thetahat)              # Number of parameters
  if(is.null(x)){
    pit = Fdist(y,thetahat,...)   # Computes the probability integral transforms
  }else{
    pit=Fdist(y,x,thetahat,...)
  }
  if(is.null(x)){
    u =t(Score(y,thetahat,...))   # Components of the score: n by p matrix
  }else{
    u =t(Score(y,x,thetahat,...)) # Components of the score: n by p matrix
  }
  Fisher = t(u)%*% u / n          # Estimate of Fisher information in 1 point.
  s=(1:m)/(m+1)                   # Grid on which to compute covariance matrix.
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
  Sigma.CvM =  m1-m2-Df%*%solve(Fisher,t(Df))
  Sigma.AD =  Sigma.CvM/sqrt(outer(s*(1-s),s*(1-s),"*"))
  J = diag(m) - matrix(1/m,m,m)
  Sigma.Watson =  J %*% Sigma.CvM %*% J
  Evals.CvM = eigen(Sigma.CvM)$values/m
  Evals.AD = eigen(Sigma.AD)$values/m
  Evals.Watson = eigen(Sigma.Watson)$values/m
  stat=gof.statistics.only(pit)
  P.CvM = imhof(stat$CvM,Evals.CvM)$Qq
  P.AD = imhof(stat$AD,Evals.AD)$Qq
  P.Watson = imhof(stat$Watson,Evals.Watson)$Qq
  #
  list(CvM=list(W2=stat$CvM,P=P.CvM),AD=list(A2=stat$AD,P=P.AD),Watson=list(U2=stat$Watson,P=P.Watson))
}
