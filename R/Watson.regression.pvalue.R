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
Watson.normal.regression.pvalue = function(w,x,neig=max(n,100),verbose=FALSE){
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
