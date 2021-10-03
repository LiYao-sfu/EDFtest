#' P-value of the EDF Statistic A^2 for a Given Distribution
#'
#' @description
#' Compute p-value of the given Anderson-Darling statistic A^2.
#'
#' @details
#'
#' @param a Anderson-Darling statistic W^2 for a given distribution.
#' @param neig Number of eigenvalues used for \code{imhof()}.
#' @param verbose Logical; if TRUE, print warning messages.
#' @param shape The shape parameter of Gamma distribution.
#'
#' @return P-value of the Anderson-Darling statistic of a uniform sample.
#' @export
#'
#' @examples
#' x1=rnorm(n=100,mean=0,sd=1)
#' asq1 = AD.normal(x1)
#' AD.normal.pvalue(asq1)
#'
#' x2=rgamma(n=100,shape=1,scale=1)
#' asq2 = AD.gamma(x)
#' AD.gamma.pvalue(asq2,1)
#'
#' x3=rlogis(n=100,location=0,scale=1)
#' asq3 = AD.logistic(x)
#' AD.logistic.pvalue(asq3)
#'
#' x4= rmutil::rlaplace(n=100,m=0,s=1)
#' asq4 = AD.laplace(x)
#' AD.laplace.pvalue(asq4)
#'
#' x5=rweibull(n=100,shape=1,scale=1)
#' asq5 = AD.weibull(x)
#' AD.weibull.pvalue(asq5)
#'
#' x6=rexp(n=100,rate=1/2)
#' asq6 = AD.exp(x)
#' AD.exp.pvalue(asq6)
AD.normal.pvalue = function(a,neig=100,verbose=FALSE){
  e = AD.normal.eigen(neig)
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
#' @rdname AD.normal.pvalue
AD.gamma.pvalue = function(a,shape,neig = 100,verbose=FALSE){
  e = AD.gamma.eigen(neig,shape=shape)
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
    if(verbose) cat("forA = ",a," and neig = ",neig,
                    " p was replaced by a lower bound on p: ",plb, "\n")
  }
  list(P=p,error=aerror)
}

#' @export
#' @rdname AD.normal.pvalue
AD.logistic.pvalue = function(a,neig=100,verbose=FALSE){
  e = AD.logistic.eigen(neig)
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
#' @rdname AD.normal.pvalue
AD.laplace.pvalue = function(a,neig=100,verbose=FALSE){
  e = AD.laplace.eigen(neig)
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
#' @rdname AD.normal.pvalue
AD.weibull.pvalue = function(a,neig=100,verbose=FALSE){
  e=AD.weibull.eigen(neig)
  plb=pchisq(a/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(a,lambda=e,epsabs = 1e-9,limit=2^7)
  options(warn=warn)
  aerror=im$abserr
  p=im$Qq
  if(p < 0 ) cat("Imhof returned a negative probability\n")
  if(p < plb){
    p=plb
    if(verbose) cat(" p =",p[i]," was replaced by a lower bound on p:",plb, "\n")
  }

  list(P=p,error=aerror)
}

#' @export
#' @rdname AD.normal.pvalue
AD.exp.pvalue = function(a,neig=100,verbose=FALSE){
  e=AD.exp.eigen(neig)
  plb=pchisq(a/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(a,lambda=e,epsabs = 1e-9,limit=2^7) #play with eps and limit
  options(warn=warn)
  aerror=im$abserr
  p=im$Qq
  if(p < 0 ) cat("Imhof returned a negative probability\n")
  if(p < plb){
    p=plb
    if(verbose) cat(" p =",p[i]," was replaced by a lower bound on p:",plb, "\n")
  }
  list(P=p,error=aerror)
}

# Helpers -----------------------------------------------------------------

AD.normal.eigen = function(n){
  M=AD.normal.covmat(n)
  e=eigen(M)$values/n
  e
}

AD.normal.covmat=function(n){
  s=1:n
  s=s/(n+1)
  t=s
  CvM.normal.covmat(n)/sqrt(outer(s*(1-s),t*(1-t)))
}

AD.gamma.eigen=function(n,shape){
  # I don't have code for this yet
  # mean = 1 -?
  s=1:n
  s=s/(n+1)
  M=AD.gamma.covmat(n,shape)
  e=eigen(M)$values/n
  e # *mean/sum(e)
}

AD.gamma.covmat=function(n,shape){
  s=1:n
  s=s/(n+1)
  CvM.gamma.covmat(n,shape=shape)/sqrt(outer(s*(1-s),s*(1-s)))
}

AD.logistic.eigen = function(n){
  mean.asq.logistic = 1-(2*pi^2-3)/(2*(pi^2+3))  # from Maple
  M=AD.logistic.covmat(n)
  e=eigen(M)$values/n
  e*mean.asq.logistic/sum(e)
}

AD.logistic.covmat=function(n){
  s=1:n
  s=s/(n+1)
  t=s
  CvM.logistic.covmat(n)/sqrt(outer(s*(1-s),t*(1-t)))
}

AD.laplace.eigen = function(n){
  mean = 1-.5351471355520514227 # from Maple
  M=AD.laplace.covmat(n)
  e=eigen(M)$values/n
  e  * mean / sum(e)
}

AD.laplace.covmat=function(n){
  s=1:n
  s=s/(n+1)
  s2=sqrt(s*(1-s))
  CvM.laplace.covmat(n)/outer(s2,s2)
}

AD.weibull.eigen = function(n){
  mean.asq.weibull = 0.3868394505  # from Maple
  M=AD.weibull.covmat(n)
  e=eigen(M)$values/n
  e*mean.asq.weibull/sum(e)
}

AD.weibull.covmat=function(n){
  s=1:n
  s=s/(n+1)
  t=s
  CvM.weibull.covmat(n)/sqrt(outer(s*(1-s),t*(1-t)))
}

AD.exp.eigen=function(n){
  mean = 0.595886194 # 3 - 2*zeta(3) from Maple
  s=1:n
  s=s/(n+1)
  M=AD.exp.covmat(n)
  e=eigen(M)$values/n
  e * mean / sum(e)
}

AD.exp.covmat=function(n){
  s=1:n
  s=s/(n+1)
  CvM.exp.covmat(n)/sqrt(outer(s*(1-s),s*(1-s)))
}

