#' P-value of EDF statistics A^2 for Exponential Distribution
#'
#' Compute p-value of the given Anderson-Darling statistic A^2
#'
#' @param a A^2 for Exponential Distribution
#' @param neig number of eigenvalues
#' @param verbose logical; if TRUE, print warning messages
#'
#' @return AD.exp.pvalue gives p-value of the Anderson-Darling statistic of a uniform sample.
#' @export
#'
#' @examples
#' x = rexp(100)
#' asq = AD.exp(x)
#' AD.exp.pvalue(asq)
AD.exp.pvalue = function(a,neig=100,verbose=FALSE){
  require("CompQuadForm")
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
