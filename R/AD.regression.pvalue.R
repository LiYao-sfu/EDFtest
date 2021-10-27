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
AD.normal.regression.pvalue = function(w,x,neig=max(n,100),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = AD.normal.regression.eigen(x,neig=neig)
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
#' @rdname AD.regression.pvalue
AD.gamma.regression.pvalue = function(w,x,theta,link="log",neig=max(n,100),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = AD.gamma.regression.eigen(x,theta=theta,link=link,neig=neig)
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
#' @rdname AD.regression.pvalue
AD.logistic.regression.pvalue = function(w,x,neig=max(n,100),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = AD.logistic.regression.eigen(x,neig=neig)
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
#' @rdname AD.regression.pvalue
AD.laplace.regression.pvalue = function(w,xneig=max(100,n),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = AD.laplace.regression.eigen(w,xneig=max(100,n))
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
#' @rdname AD.regression.pvalue
AD.weibull.regression.pvalue = function(w,x,neig=max(n,100),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = AD.weibull.regression.eigen(x,neig=neig)
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
#' @rdname AD.regression.pvalue
AD.extremevalue.regression.pvalue = function(w,x,neig=max(n,100),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = AD.weibull.regression.eigen(x,neig=neig)
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
#' @rdname AD.regression.pvalue
AD.exp.regression.pvalue = function(w,x,theta,link="log",neig=max(n,100),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = AD.exp.regression.eigen(x,theta=theta,link=link,neig=neig)
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


