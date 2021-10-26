#' P-value of Cramer-von Mises statistic for regression
#'
#' @description
#'
#' @param w Cramér-von Mises statistic \eqn{W^2} with a given distribution.
#' @param x explanatory variables
#' @param neig
#' @param verbose
#'
#' @return
#'
#' @seealso
#'
#' @name CvM.regression.pvalue
#' @examples
#'
NULL

#' @export
#' @rdname CvM.regression.pvalue
CvM.normal.regression.pvalue = function(w,x,neig=max(n,100),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = CvM.normal.regression.eigen(x,neig=neig)
  plb=pchisq(w/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(w,lambda=e,epsabs = 1e-9,limit=2^7)
  options(warn=warn)
  aerror=im$abserr
  p=im$Qq
  if(p<0&&verbose)cat("for W = ",w," using = ",neig,
                      " eigenvalues, imhof returned a negative probability\n")
  if(p<plb){
    p=plb
    if(verbose) cat("for W = ",w," using = ",neig,
                    " eigenvalues, p was replaced by a lower bound on p: ",plb, "\n")
  }
  list(P=p,error=aerror)
}



#' @export
#' @rdname CvM.regression.pvalue
CvM.gamma.regression.pvalue = function(w,x,neig=max(n,100),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = CvM.gamma.regression.eigen(x,neig=neig)
  plb=pchisq(w/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(w,lambda=e,epsabs = 1e-9,limit=2^7)
  options(warn=warn)
  aerror=im$abserr
  p=im$Qq
  if(p<0&&verbose)cat("for W = ",w," using = ",neig,
                      " eigenvalues, imhof returned a negative probability\n")
  if(p<plb){
    p=plb
    if(verbose) cat("for W = ",w," using = ",neig,
                    " eigenvalues, p was replaced by a lower bound on p: ",plb, "\n")
  }
  list(P=p,error=aerror)
}


#' @export
#' @rdname CvM.regression.pvalue
CvM.logistic.regression.pvalue = function(w,x,neig=max(n,100),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = CvM.logistic.regression.eigen(x,neig=neig)
  plb=pchisq(w/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(w,lambda=e,epsabs = 1e-9,limit=2^7)
  options(warn=warn)
  aerror=im$abserr
  p=im$Qq
  if(p<0&&verbose)cat("for W = ",w," using = ",neig,
                      " eigenvalues, imhof returned a negative probability\n")
  if(p<plb){
    p=plb
    if(verbose) cat("for W = ",w," using = ",neig,
                    " eigenvalues, p was replaced by a lower bound on p: ",plb, "\n")
  }
  list(P=p,error=aerror)
}


#' @export
#' @rdname CvM.regression.pvalue
CvM.laplace.regression.pvalue = function(w,xneig=max(100,n),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = CvM.laplace.regression.eigen(w,xneig=max(100,n))
  plb=pchisq(w/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(w,lambda=e,epsabs = 1e-9,limit=2^7)
  options(warn=warn)
  aerror=im$abserr
  p=im$Qq
  if(p<0&&verbose)cat("for W = ",w," and neig = ",neig,
                      " imhof returned a negative probability\n")
  if(p<plb){
    p=plb
    if(verbose) cat("for W = ",w," and neig = ",neig,
                    " p was replaced by a lower bound on p: ",plb, "\n")
  }
  list(P=p,error=aerror)
}


#' @export
#' @rdname CvM.regression.pvalue
CvM.weibull.regression.pvalue = function(w,x,neig=max(n,100),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = CvM.weibull.regression.eigen(x,neig=neig)
  plb=pchisq(w/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(w,lambda=e,epsabs = 1e-9,limit=2^7)
  options(warn=warn)
  aerror=im$abserr
  p=im$Qq
  if(p<0&&verbose)cat("for W = ",w," using = ",neig,
                      " eigenvalues, imhof returned a negative probability\n")
  if(p<plb){
    p=plb
    if(verbose) cat("for W = ",w," using = ",neig,
                    " eigenvalues, p was replaced by a lower bound on p: ",plb, "\n")
  }
  list(P=p,error=aerror)
}


#' @export
#' @rdname CvM.regression.pvalue
CvM.extremevalue.regression.pvalue = function(w,x,neig=max(n,100),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = CvM.weibull.regression.eigen(x,neig=neig)
  plb=pchisq(w/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(w,lambda=e,epsabs = 1e-9,limit=2^7)
  options(warn=warn)
  aerror=im$abserr
  p=im$Qq
  if(p<0&&verbose)cat("for W = ",w," using = ",neig,
                      " eigenvalues, imhof returned a negative probability\n")
  if(p<plb){
    p=plb
    if(verbose) cat("for W = ",w," using = ",neig,
                    " eigenvalues, p was replaced by a lower bound on p: ",plb, "\n")
  }
  list(P=p,error=aerror)
}


#' @export
#' @rdname CvM.regression.pvalue
CvM.exp.regression.pvalue = function(w,x,neig=max(n,100),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = CvM.exp.regression.eigen(x,neig=neig)
  plb=pchisq(w/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(w,lambda=e,epsabs = 1e-9,limit=2^7)
  options(warn=warn)
  aerror=im$abserr
  p=im$Qq
  if(p<0&&verbose)cat("for W = ",w," using = ",neig,
                      " eigenvalues, imhof returned a negative probability\n")
  if(p<plb){
    p=plb
    if(verbose) cat("for W = ",w," using = ",neig,
                    " eigenvalues, p was replaced by a lower bound on p: ",plb, "\n")
  }
  list(P=p,error=aerror)
}

# Helpers -----------------------------------------------------------------


CvM.normal.regression.eigen = function(x,neig=max(n,100)){
  p=dim(x)[2]
  n=dim(x)[1]
  mean.wsq.normal=1/6 -7*sqrt(3)/(36*pi)
  M=CvM.normal.regression.covmat(x,neig=neig)
  e=eigen(M)$values/neig
  e*mean.wsq.normal/sum(e)
}


CvM.normal.regression.covmat=function(x,m=max(n,100)){
  p=dim(x)[2]
  n = dim(x)[1]
  D = t(x)%*%x/n
  Fisher.normal = matrix(0,nrow=p+1,ncol=p+1)
  Fisher.normal[1:p,1:p]= D
  Fisher.normal[p+1,p+1]=2
  s=1:m
  s=s/(m+1)
  M1=outer(s,s,pmin)-outer(s,s)
  x = qnorm(s)
  G1 = dnorm(x)
  G2 = -x*G1
  G1 = x*rep(G1,p)
  M2 = cbind(G1,G2)
  M1-M2%*%solve(Fisher.normal,t(M2))
}


CvM.gamma.regression.covmat=function(x,theta,link="log",neig=max(n,100)){
  fisher.information.gamma=function(x,theta){
    #
    # returns the estimated Fisher Information per point
    # for a gamma regression model in which the log mean is predicted
    # linearly from a matrix of covariates x
    # Normally x will contain an intercept term
    #
    pp=length(theta)
    p=pp-1
    shape.hat = theta[pp]
    coefs = theta[-pp]
    eta = x %*% coefs
    if( link == "log") {
      invlink = exp
      linkder = exp
    }
    if( link == "inverse" ) {
      invlink = function(w) 1/w
      linkder = function(w) -1/w^2
    }
    if( link == "identity" ) {
      invlink = function(w) w
      lindker = function(w) 1
    }
    mu = invlink(eta)
    deriv = linkder(eta)
    W = x * rep(deriv,p)
    M = t(W)%*%W/n  # should be p by p
    print(dim(M))
    FI=matrix(0,nrow=pp,ncol=pp)
    FI[pp,pp]=trigamma(shape.hat)-1/shape.hat
    FI[1:p,1:p]=(shape.hat/mu^2) * M
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



CvM.gamma.regression.eigen = function(x,theta,link="log",neig=max(n,100)){
  p=dim(x)[2]
  n=dim(x)[1]
  M=CvM.gamma.regression.covmat(x,theta=theta,link=link,neig=neig)
  e=eigen(M)$values/neig
  e # *mean.wsq.gamma/sum(e) but this correction is not available yet
}


CvM.logistic.regression.eigen = function(x,neig=max(100,n)){
  p = dim(x)[2]
  n = dim(x)[1]
  mean.wsq.logistic= 1/6 -(4*pi^2-9)/(20*(pi^2+3))  # from Maple
  M=CvM.logistic.regression.covmat(x,neig=neig)
  e=eigen(M)$values/neig
  e*mean.wsq.logistic/sum(e)
}


CvM.logistic.regression.covmat=function(x,neig=max(100,n)){
  p=dim(x)[2]
  n = dim(x)[1]
  D = t(x)%*%x/n
  Fisher = matrix(0,nrow=p+1,ncol=p+1)
  Fisher[1:p,1:p] = D/3
  Fisher[p+1,p+1] = (pi^2+3)/9
  s=1:neig
  s=s/(neig+1)
  M1=outer(s,s,pmin)-outer(s,s)
  G1 = s*(1-s)
  G2 = G1*(log(s/(1-s)))
  G1 = x*rep(G1,p)
  M2 = cbind(G1,G2)
  M1-M2%*%solve(Fisher,t(M2))
}


CvM.laplace.regression.eigen  = function(x,neig = max(n,100)){
  p = dim(x)[2]
  n = dim(x)[1]
  mean = 1/6 -1/12-1/54
  M=CvM.laplace.regression.covmat(x,neig = neig)
  e=eigen(M)$values/neig
  e * mean / sum(e)
}


CvM.laplace.regression.covmat=function(x,neig = max(n,100)){
  p=dim(x)[2]
  n = dim(x)[1]
  D = t(x)%*%x/n
  Fisher = matrix(0,nrow=p+1,ncol=p+1)
  Fisher[1:p,1:p] = D
  Fisher[p+1,p+1] = 1
  s=1:neig
  s=s/(neig+1)
  M1=outer(s,s,pmin)-outer(s,s)
  G1 = -s
  G1[s>0.5]=(s-1)[s>0.5]
  G2 = -s*log(2*s)
  G2[s>0.5] = ((1-s)*log(2*(1-s)))[s>0.5]
  G1 = x*rep(G1,p)
  M2 = cbind(G1,G2)
  M1-M2%*%solve(Fisher.laplace,t(M2))
}


CvM.weibull.regression.eigen = function(x,neig=max(100,n)){
  p = dim(x)[2]
  n = dim(x)[1]
  # mean.wsq.weibull= 1/6 -(4*pi^2-9)/(20*(pi^2+3))  # from Maple
  M=CvM.weibull.regression.covmat(x,neig=neig)
  e=eigen(M)$values/neig
  e #*mean.wsq.logistic/sum(e)
}


CvM.weibull.regression.covmat=function(x,neig=max(100,n)){
  p=dim(x)[2]
  n = dim(x)[1]
  D = t(x)%*%x/n
  Fisher.weibull = matrix(c(pi^2/6+(1+digamma(1))^2,-(1+digamma(1)),-(1+digamma(1)),1),nrow=2)
  Fisher = matrix(0,nrow=p+1,ncol=p+1)
  Fisher[1:p,1:p] = pi^2/6+(1+digamma(1))^2*D
  Fisher[1:p,p+1] = -(1+digamma(1))*apply(x,2,mean)
  Fisher[p+1,p+1] = 1
  s=1:neig
  s=s/(neig+1)
  M1=outer(s,s,pmin)-outer(s,s)
  G1 = s*(1-s)
  G2 = G1*(log(s/(1-s)))
  G1 = x*rep(G1,p)
  M2 = cbind(G1,G2)
  M1-M2%*%solve(Fisher,t(M2))
}
