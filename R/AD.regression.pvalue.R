#' P-value of Anderson-Darling statistic for regression
#'
#' @description
#'
#' @param w Anderson-Darling statistic \eqn{A^2} with a given distribution.
#' @param x explanatory variables
#' @param neig
#' @param verbose
#'
#' @return
#'
#' @seealso
#'
#' @name AD.regression.pvalue
#' @examples
#'
NULL

#' @export
#' @rdname AD.regression.pvalue
AD.normal.regression.pvalue = function(a,x,neig=max(n,400),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = AD.normal.regression.eigen(x,neig=neig)
  plb=pchisq(a/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(a,lambda=e,epsabs = 1e-9,limit=2^7)
  options(warn=warn)
  aerror=im$abserr
  p=im$Qq
  if(p<0&&verbose)cat("for A = ",a," using = ",neig,
                      " eigenvalues, imhof returned a negative probability\n")
  if(p<plb){
    p=plb
    if(verbose) cat("for A = ",a," using = ",neig,
                    " eigenvalues, p was replaced by a lower bound on p: ",plb, "\n")
  }
  list(P=p,error=aerror)
}


#' @export
#' @rdname AD.regression.pvalue
AD.gamma.regression.pvalue = function(a,x,theta,link="log",neig=max(n,400),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = AD.gamma.regression.eigen(x,theta=theta,link=link,neig=neig)
  plb=pchisq(a/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(a,lambda=e,epsabs = 1e-9,limit=2^7)
  options(warn=warn)
  aerror=im$abserr
  p=im$Qq
  if(p<0&&verbose)cat("for A = ",a," using = ",neig,
                      " eigenvalues, imhof returned a negative probability\n")
  if(p<plb){
    p=plb
    if(verbose) cat("for A = ",a," using = ",neig,
                    " eigenvalues, p was replaced by a lower bound on p: ",plb, "\n")
  }
  list(P=p,error=aerror)
}


#' @export
#' @rdname AD.regression.pvalue
AD.logistic.regression.pvalue = function(a,x,neig=max(n,400),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = AD.logistic.regression.eigen(x,neig=neig)
  plb=pchisq(a/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(a,lambda=e,epsabs = 1e-9,limit=2^7)
  options(warn=warn)
  aerror=im$abserr
  p=im$Qq
  if(p<0&&verbose)cat("for A = ",a," using = ",neig,
                      " eigenvalues, imhof returned a negative probability\n")
  if(p<plb){
    p=plb
    if(verbose) cat("for A = ",a," using = ",neig,
                    " eigenvalues, p was replaced by a lower bound on p: ",plb, "\n")
  }
  list(P=p,error=aerror)
}


#' @export
#' @rdname AD.regression.pvalue
AD.laplace.regression.pvalue = function(a,x,neig=max(400,n),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = AD.laplace.regression.eigen(x,neig=neig)
  plb=pchisq(a/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(a,lambda=e,epsabs = 1e-9,limit=2^7)
  options(warn=warn)
  aerror=im$abserr
  p=im$Qq
  if(p<0&&verbose)cat("for A = ",a," and neig = ",neig,
                      " imhof returned a negative probability\n")
  if(p<plb){
    p=plb
    if(verbose) cat("for A = ",a," and neig = ",neig,
                    " p was replaced by a lower bound on p: ",plb, "\n")
  }
  list(P=p,error=aerror)
}


#' @export
#' @rdname AD.regression.pvalue
AD.weibull.regression.pvalue = function(a,x,neig=max(n,400),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = AD.weibull.regression.eigen(x,neig=neig)
  plb=pchisq(a/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(a,lambda=e,epsabs = 1e-9,limit=2^7)
  options(warn=warn)
  aerror=im$abserr
  p=im$Qq
  if(p<0&&verbose)cat("for A = ",a," using = ",neig,
                      " eigenvalues, imhof returned a negative probability\n")
  if(p<plb){
    p=plb
    if(verbose) cat("for A = ",a," using = ",neig,
                    " eigenvalues, p was replaced by a lower bound on p: ",plb, "\n")
  }
  list(P=p,error=aerror)
}


#' @export
#' @rdname AD.regression.pvalue
AD.extremevalue.regression.pvalue = function(a,x,neig=max(n,400),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = AD.weibull.regression.eigen(x,neig=neig)
  plb=pchisq(a/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(a,lambda=e,epsabs = 1e-9,limit=2^7)
  options(warn=warn)
  aerror=im$abserr
  p=im$Qq
  if(p<0&&verbose)cat("for A = ",a," using = ",neig,
                      " eigenvalues, imhof returned a negative probability\n")
  if(p<plb){
    p=plb
    if(verbose) cat("for A = ",a," using = ",neig,
                    " eigenvalues, p was replaced by a lower bound on p: ",plb, "\n")
  }
  list(P=p,error=aerror)
}


#' @export
#' @rdname AD.regression.pvalue
AD.exp.regression.pvalue = function(a,x,theta,link="log",neig=max(n,400),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = AD.exp.regression.eigen(x,theta=theta,link=link,neig=neig)
  plb=pchisq(a/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(a,lambda=e,epsabs = 1e-9,limit=2^7)
  options(warn=warn)
  aerror=im$abserr
  p=im$Qq
  if(p<0&&verbose)cat("for A = ",a," using = ",neig,
                      " eigenvalues, imhof returned a negative probability\n")
  if(p<plb){
    p=plb
    if(verbose) cat("for A = ",a," using = ",neig,
                    " eigenvalues, p was replaced by a lower bound on p: ",plb, "\n")
  }
  list(P=p,error=aerror)
}

# Helpers -----------------------------------------------------------------


AD.normal.regression.eigen = function(x,neig){
  p=dim(x)[2]
  n=dim(x)[1]
  mean.wsq.normal=1/6 -7*sqrt(3)/(36*pi)
  M=AD.normal.regression.covmat(x,neig=neig)
  e=eigen(M)$values/neig
  e*mean.wsq.normal/sum(e)
}

# CvM.normal.regression.covmat=function(x,neig){
#   p=dim(x)[2]
#   n = dim(x)[1]
#   D = t(x)%*%x/n
#   Fisher.normal = matrix(0,nrow=p+1,ncol=p+1)
#   Fisher.normal[1:p,1:p]= D
#   Fisher.normal[p+1,p+1]=2
#   s=1:neig
#   s=s/(neig+1)
#   M1=outer(s,s,pmin)-outer(s,s)
#   xx = qnorm(s)
#   G1 = dnorm(xx)
#   G2 = -xx*G1
#   Del = apply(x,2,mean) # Del should now be a p vector
#   G1 = outer(Del, G1) # G1 should now be p by neig
#   M2 = rbind(G1,G2) # M2 should now be p+1 by neig
#   M1 - t(M2) %*% solve(Fisher.normal,M2)
# }

AD.normal.regression.covmat=function(x,neig){
  p=dim(x)[2]
  n = dim(x)[1]
  D = t(x)%*%x/n
  Fisher.normal = matrix(0,nrow=p+1,ncol=p+1)
  Fisher.normal[1:p,1:p]= D
  Fisher.normal[p+1,p+1]=2
  s=1:neig
  s=s/(neig+1)
  w = 1/sqrt(s*(1-s))
  M1=outer(s,s,pmin)-outer(s,s)
  xx = qnorm(s)
  G1 = dnorm(xx)
  G2 = -xx*G1
  Del = apply(x,2,mean) # Del should now be a p vector
  G1 = outer(Del, G1)
  #G1 = xx*rep(G1,p)
  M2 = rbind(G1,G2)
  (M1-M2%*%solve(Fisher.normal,t(M2)))/outer(w,w)
}

AD.gamma.regression.eigen = function(x,theta,neig,link="log"){
  p=dim(x)[2]
  n=dim(x)[1]
  M=AD.gamma.regression.covmat(x,theta=theta,link=link,neig=neig)
  e=eigen(M)$values/neig
  e # *mean.wsq.gamma/sum(e) but this correction is not available yet
}


AD.gamma.regression.covmat=function(x,theta,neig,link="log"){
  M.CvM = CvM.gamma.regression.covmat(x,theta,link=link,neig=neig)
  s=1:neig
  s=s/(neig+1)
  w = 1/sqrt(s*(1-s))
  return(M.CvM/outer(w,w))
}

#Delete?
AD.gamma.regression.covmat.old=function(x,theta,neig,link="log"){
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
    W = x * rep(deriv/mu,p)
    FI=matrix(0,nrow=pp,ncol=pp)
    FI[pp,pp]=trigamma(shape.hat)-1/shape.hat
    FI[1:p,1:p] = t(W)%*%W/n  # should be p by p
    FI
  }
  g = gamma(shape)
  dg = digamma(shape)
  FI = fisher.information.gamma(shape.hat = shape)
  s=1:n
  s=s/(n+1)
  w = 1/sqrt(s*(1-s))
  w = 1/sqrt(s*(1-s))
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
  (M1-M2%*%solve(FI,t(M2)))/outer(w,w)
}


AD.exp.regression.eigen = function(x,theta,link="log",neig=max(n,400)){
  p=dim(x)[2]
  n=dim(x)[1]
  M=AD.exp.regression.covmat(x,theta=theta,link=link,neig=neig)
  e=eigen(M)$values/neig
  e # *mean.wsq.exp/sum(e) but this correction is not available yet
}


AD.exp.regression.covmat=function(x,theta,neig,link="log"){
  M.CvM = CvM.exp.regression.covmat(x,theta,link=link,neig=neig)
  s=1:neig
  s=s/(neig+1)
  w = 1/sqrt(s*(1-s))
  return(M.CvM/outer(w,w))
}

# AD.exp.regression.covmat=function(x,theta,link="log",neig=max(n,400)){
#     #
#     # returns the estimated Fisher Information per point
#     # for an exponential regression model in which the log mean is predicted
#     # linearly from a matrix of covariates x
#     # Normally x will contain an intercept term
#     #
#     pp=length(theta)
#     eta = x %*% theta
#     if( link == "log") {
#       invlink = exp
#       linkder = exp
#     }
#     if( link == "inverse" ) {
#       invlink = function(w) 1/w
#       linkder = function(w) -1/w^2
#     }
#     if( link == "identity" ) {
#       invlink = function(w) w
#       lindker = function(w) 1
#     }
#     mu = invlink(eta)
#     deriv = linkder(eta)
#     W = x * rep(deriv/mu,p)
#     FI =  t(W)%*%W/n  # should be p by p
#   s = 1:n
#   s = s/(n+1)
#   M1 = outer(s,s,pmin)-outer(s,s)
#   Q = qgamma(s,shape=shape)
#   D = dgamma(Q,shape=shape)
#   M2 = - Q * D
#   (M1-M2%*%solve(FI,t(M2)))/outer(w,w)
# }

AD.logistic.regression.eigen = function(x,neig){
  p = dim(x)[2]
  n = dim(x)[1]
  mean.wsq.logistic= 1/6 -(4*pi^2-9)/(20*(pi^2+3))  # from Maple
  M=AD.logistic.regression.covmat(x,neig=neig)
  e=eigen(M)$values/neig
  e*mean.wsq.logistic/sum(e)
}


AD.logistic.regression.covmat=function(x,neig){
  M.CvM = CvM.logistic.regression.covmat(x,neig=neig)
  s=1:neig
  s=s/(neig+1)
  w = 1/sqrt(s*(1-s))
  return(M.CvM/outer(w,w))
}



AD.logistic.regression.covmat.old=function(x,neig){
  p=dim(x)[2]
  n = dim(x)[1]
  D = t(x)%*%x/n
  Fisher = matrix(0,nrow=p+1,ncol=p+1)
  Fisher[1:p,1:p] = D/3
  Fisher[p+1,p+1] = (pi^2+3)/9
  s=1:neig
  s=s/(neig+1)
    w = 1/sqrt(s*(1-s))
  M1=outer(s,s,pmin)-outer(s,s)
  G1 = s*(1-s)
  G2 = G1*(log(s/(1-s)))
  G1 = x*rep(G1,p)
  M2 = cbind(G1,G2)
  (M1-M2%*%solve(FI,t(M2)))/outer(w,w)
}


AD.laplace.regression.eigen  = function(x,neig){
  p = dim(x)[2]
  n = dim(x)[1]
#  mean = 1/6 -1/12-1/54
  M=AD.laplace.regression.covmat(x,neig = neig)
  e=eigen(M)$values/neig
  e  # * mean / sum(e)
}

AD.laplace.regression.covmat=function(x,neig){
  M.CvM = CvM.laplace.regression.covmat(x,neig=neig)
  s=1:neig
  s=s/(neig+1)
  w = 1/sqrt(s*(1-s))
  return(M.CvM/outer(w,w))
}

AD.laplace.regression.covmat.old=function(x,neig){
  p=dim(x)[2]
  n = dim(x)[1]
  D = t(x)%*%x/n
  Fisher = matrix(0,nrow=p+1,ncol=p+1)
  Fisher[1:p,1:p] = D
  Fisher[p+1,p+1] = 1
  s=1:neig
  s=s/(neig+1)
  w = 1/sqrt(s*(1-s))
  M1=outer(s,s,pmin)-outer(s,s)
  G1 = -s
  G1[s>0.5]=(s-1)[s>0.5]
  G2 = -s*log(2*s)
  G2[s>0.5] = ((1-s)*log(2*(1-s)))[s>0.5]
  G1 = x*rep(G1,p)
  M2 = cbind(G1,G2)
  (M1-M2%*%solve(FI,t(M2)))/outer(w,w)
}


AD.weibull.regression.eigen = function(x,neig){
  p = dim(x)[2]
  n = dim(x)[1]
  # mean.wsq.weibull= 1/6 -(4*pi^2-9)/(20*(pi^2+3))  # from Maple
  M=AD.weibull.regression.covmat(x,neig=neig)
  e=eigen(M)$values/neig
  e #*mean.wsq.logistic/sum(e)
}


AD.weibull.regression.covmat=function(x,neig){
  M.CvM = CvM.weibull.regression.covmat(x,neig=neig)
  s=1:neig
  s=s/(neig+1)
  w = 1/sqrt(s*(1-s))
  return(M.CvM/outer(w,w))
}


AD.weibull.regression.covmat.old=function(x,neig=max(400,n)){
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
    w = 1/sqrt(s*(1-s))
  M1=outer(s,s,pmin)-outer(s,s)
  G1 = s*(1-s)
  G2 = G1*(log(s/(1-s)))
  G1 = x*rep(G1,p)
  M2 = cbind(G1,G2)
  (M1-M2%*%solve(FI,t(M2)))/outer(w,w)
}

