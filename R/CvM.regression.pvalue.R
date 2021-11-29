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
CvM.normal.regression.pvalue = function(w,x,neig=max(n,400),verbose=FALSE){
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
                    " eigenvalues, p-value was replaced by a lower bound: ",plb, "\n")
  }
  list(P=p,error=aerror)
}


#' @export
#' @rdname CvM.regression.pvalue
CvM.gamma.regression.pvalue = function(w,x,theta,neig=max(n,400),link="log",verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = CvM.gamma.regression.eigen(x,theta=theta,link=link,neig=neig)
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
                    " eigenvalues, p-value was replaced by a lower bound: ",plb, "\n")
  }
  list(P=p,error=aerror)
}


#' @export
#' @rdname CvM.regression.pvalue
CvM.logistic.regression.pvalue = function(w,x,neig=max(n,400),verbose=FALSE){
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
                    " eigenvalues, p-value was replaced by a lower bound: ",plb, "\n")
  }
  list(P=p,error=aerror)
}


#' @export
#' @rdname CvM.regression.pvalue
CvM.laplace.regression.pvalue = function(w,x,neig=max(400,n),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = CvM.laplace.regression.eigen(x,neig=neig)
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
                    " p-value was replaced by a lower bound: ",plb, "\n")
  }
  list(P=p,error=aerror)
}


#' @export
#' @rdname CvM.regression.pvalue
CvM.weibull.regression.pvalue = function(w,x,neig=max(n,400),verbose=FALSE){
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
                    " eigenvalues, p-value was replaced by a lower bound: ",plb, "\n")
  }
  list(P=p,error=aerror)
}


#' @export
#' @rdname CvM.regression.pvalue
CvM.extremevalue.regression.pvalue = function(w,x,neig=max(n,400),verbose=FALSE){
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
                    " eigenvalues, p-value was replaced by a lower bound: ",plb, "\n")
  }
  list(P=p,error=aerror)
}


#' @export
#' @rdname CvM.regression.pvalue
CvM.exp.regression.pvalue = function(w,x,theta,neig=max(n,400),link="log",verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = CvM.exp.regression.eigen(x,theta=theta,link=link,neig=neig)
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
                    " eigenvalues, p-value was replaced by a lower bound: ",plb, "\n")
  }
  list(P=p,error=aerror)
}

# Helpers -----------------------------------------------------------------


CvM.normal.regression.eigen = function(x,neig){
  mean.wsq.normal=1/6 -7*sqrt(3)/(36*pi)
  M=CvM.normal.regression.covmat(x,neig=neig)
  e=eigen(M)$values/neig
  e*mean.wsq.normal/sum(e)
}


CvM.normal.regression.covmat=function(x,neig){
  p=dim(x)[2]
  n = dim(x)[1]
  D = t(x)%*%x/n
  Fisher.normal = matrix(0,nrow=p+1,ncol=p+1)
  Fisher.normal[1:p,1:p]= D
  Fisher.normal[p+1,p+1]=2
  s=1:neig
  s=s/(neig+1)
  M1=outer(s,s,pmin)-outer(s,s)
  xx = qnorm(s)
  G1 = dnorm(xx)
  G2 = -xx*G1
  Del = apply(x,2,mean) # Del should now be a p vector
  G1 = outer(Del, G1) # G1 should now be p by neig
  M2 = rbind(G1,G2) # M2 should now be p+1 by neig
  M1 - t(M2) %*% solve(Fisher.normal,M2)
}


CvM.gamma.regression.covmat=function(x,theta,neig,link="log"){
    if( link == "log") {
      invlink = exp
      linkder = exp
      weight = function(w) 1
    }
    if( link == "inverse" ) {
      invlink = function(w) 1/w
      linkder = function(w) -1/w^2
      weight = function(w) -1/w
    }
    if( link == "identity" ) {
      invlink = function(w) w
      lindker = function(w) 1
      weight = function(w) 1/w
    }
    pp=length(theta)
    p=pp-1
    shape = theta[pp]
    coefs = theta[-pp]
    eta = x %*% coefs

    mu = invlink(eta)
    wt = weight(eta)
    W = x * rep(wt,p)
    M = t(W)%*%W/n  # should be p by p
  #####
    #
    # the estimated Fisher Information per point
    # for a gamma regression model in which the log mean is predicted
    # linearly from a matrix of covariates x
    # Normally x will contain an intercept term
    #
    Fisher=matrix(0,nrow=pp,ncol=pp)
    Fisher[pp,pp]=trigamma(shape)-1/shape
    Fisher[1:p,1:p] = shape * M
    #  
  #####
  
  s=1:neig
  s=s/(neig+1)
  M1=outer(s,s,pmin)-outer(s,s)

  Q = qgamma(s,shape=shape)
  D = dgamma(Q,shape=shape)
  G1base = - Q * D
  Del = apply(W,2,mean) # Del should now be a p vector
  G1 = outer(Del, G1base) # G1 should now be p by neig

  g = gamma(shape)
  dg = digamma(shape)
  G2 = s*0
  g2.integrand = function(x,shape){log(x)*x^(shape-1)*exp(-x)}
  for(i in 1:neig){
    G2[i] = integrate(g2.integrand,0,Q[i],shape=shape)$value/g -s[i]*dg
  }
  G2=G2+G1base/shape
  M2 = rbind(G1,G2) # M2 should be p+1 by neig
  M1 - t(M2) %*% solve(Fisher,M2)
}



CvM.gamma.regression.eigen = function(x,theta,neig,link="log"){
  M=CvM.gamma.regression.covmat(x,theta=theta,link=link,neig=neig)
  e=eigen(M)$values/neig
  e # *mean.wsq.gamma/sum(e) but this correction is not available yet
}


CvM.exp.regression.covmat=function(x,theta,neig,link="log"){
    #
    # for an exponential regression model in which the log mean is predicted
    # linearly from a matrix of covariates x
    # Normally x will contain an intercept term
    # Permits use of two other links: identity and inverse
    #
  p=length(theta)
  eta = x %*% theta
    if( link == "log") {
   #  invlink = exp
   #  linkder = exp
      weight = function(w) 1
    }
    if( link == "inverse" ) {
   #  invlink = function(w) 1/w
   #  linkder = function(w) -1/w^2
      weight = function(w) -1/w
    }
    if( link == "identity" ) {
   #  invlink = function(w) w
   #  lindker = function(w) 1
      weight = function(w) 1/w
    }
  wt = weight(eta)
  W = x * rep(wt,p)
  Fisher = t(W)%*%W/n  # should be p by p
  s = 1:neig
  s = s/(neig+1)
  M1 = outer(s,s,pmin)-outer(s,s)
  G1 = -log(1-s)*(1-s)
  Del = apply(W,2,mean) # Del should now be a p vector
  G1 = outer(Del, G1) # G1 should now be p by neig
  M1 - t(G1) %*% solve(Fisher,G1)
}



CvM.exp.regression.eigen = function(x,theta,neig,link="log"){
  M=CvM.exp.regression.covmat(x,theta=theta,link=link,neig=neig)
  e=eigen(M)$values/neig
  e # *mean.wsq.exp/sum(e) but this correction is not available yet
}


CvM.logistic.regression.eigen = function(x,neig){
  mean.wsq.logistic= 1/6 -(4*pi^2-9)/(20*(pi^2+3))  # from Maple
  M=CvM.logistic.regression.covmat(x,neig=neig)
  e=eigen(M)$values/neig
  e*mean.wsq.logistic/sum(e)
}


CvM.logistic.regression.covmat=function(x,neig){
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
  Del = apply(x,2,mean) # Del should now have length p
  G1 = outer(Del,G1) # G1 should now have dimension p by neig
  M2 = rbind(G1,G2) # M2 should now have dimension p+1 by neig
  M1-t(M2)%*%solve(Fisher,M2)
}


CvM.laplace.regression.eigen  = function(x,neig){
  mean = 1/6 -1/12-1/54
  M=CvM.laplace.regression.covmat(x,neig = neig)
  e=eigen(M)$values/neig
  e * mean / sum(e)
}


CvM.laplace.regression.covmat=function(x,neig){
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
  Del = apply(x,2,mean) # Del should now have length p
  G1 = outer(Del,G1) # G1 should now have dimension p by neig
  M2 = rbind(G1,G2) # M2 should now have dimension p+1 by neig
  M1 - t(M2) %*% solve(Fisher,M2)
}


CvM.weibull.regression.eigen = function(x,neig){
  mean.wsq.weibull= 1/54 -4*(log(3)^2-log(3)-1)/(9*pi^2)  # from Maple
  M=CvM.weibull.regression.covmat(x,neig=neig)
  e=eigen(M)$values/neig
  e * mean.wsq./sum(e)
}


CvM.weibull.regression.covmat=function(x,neig){
  p=dim(x)[2]
  n = dim(x)[1]
  D = t(x)%*%x/n
  Fisher.weibull = matrix(c(pi^2/6+(1+digamma(1))^2,-(1+digamma(1)),-(1+digamma(1)),1),nrow=2)
  Fisher = matrix(0,nrow=p+1,ncol=p+1)
  Fisher[1:p,1:p] = D
  Fisher[1:p,p+1] = -(1+digamma(1))*apply(x,2,mean)
  Fisher[p+1,1:p] = -(1+digamma(1))*apply(x,2,mean)
  Fisher[p+1,p+1] = pi^2/6+(1+digamma(1))^2
  s=1:neig
  s=s/(neig+1)
  M1=outer(s,s,pmin)-outer(s,s)
  G1 = log(1-s)*(1-s)
  G2 = -log(1-s)*log(-log(1-s))*(1-s)
  Del = apply(x,2,mean) # Del should now have length p
  G1 = outer(Del,G1) # G1 should now have dimension p by neig
  M2 = rbind(G1,G2) # M2 should now have dimension p+1 by neig
  M1 - t(M2) %*% solve(Fisher,M2)
}
