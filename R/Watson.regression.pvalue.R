#' P-value of Watson statistic for regression
#'
#' @description
#'
#' @param w Watson statistic \eqn{U^2} with a given distribution.
#' @param x explanatory variables
#' @param neig
#' @param verbose
#'
#' @return
#'
#' @seealso
#'
#' @name Watson.regression.pvalue
#' @examples
#'
NULL

#' @export
#' @rdname Watson.regression.pvalue
Watson.normal.regression.pvalue = function(w,x,neig=max(n,400),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = Watson.normal.regression.eigen(x,neig=neig)
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
#' @rdname Watson.regression.pvalue
Watson.gamma.regression.pvalue = function(w,x,theta,link="log",neig=max(n,400),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = Watson.gamma.regression.eigen(x,theta=theta,link=link,neig=neig)
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
#' @rdname Watson.regression.pvalue
Watson.logistic.regression.pvalue = function(w,x,neig=max(n,400),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = Watson.logistic.regression.eigen(x,neig=neig)
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
#' @rdname Watson.regression.pvalue
Watson.laplace.regression.pvalue = function(w,xneig=max(400,n),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = Watson.laplace.regression.eigen(w,xneig=max(100,n))
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
#' @rdname Watson.regression.pvalue
Watson.weibull.regression.pvalue = function(w,x,neig=max(n,400),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = Watson.weibull.regression.eigen(x,neig=neig)
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
#' @rdname Watson.regression.pvalue
Watson.extremevalue.regression.pvalue = function(w,x,neig=max(n,400),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = Watson.weibull.regression.eigen(x,neig=neig)
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
#' @rdname Watson.regression.pvalue
Watson.exp.regression.pvalue = function(w,x,theta,link="log",neig=max(n,400),verbose=FALSE){
  p=dim(x)[2]
  n=dim(x)[1]
  e = Watson.exp.regression.eigen(x,theta=theta,link=link,neig=neig)
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


