#' P-value of EDF statistics A^2 for Weibull Distribution
#'
#' Compute p-value of the given Anderson-Darling statistic A^2
#'
#' @param a A^2 for Weibull Distribution
#' @param neig number of eigenvalues
#' @param verbose logical; if TRUE, print warning messages
#'
#' @return AD.weibull.pvalue gives p-value of the Anderson-Darling statistic of a uniform sample.
#' @export
#'
#' @examples
#' x = rweibull(1000,1)
#' asq = AD.weibull(x)
#' AD.weibull.pvalue(asq)
AD.weibull.pvalue = function(a,neig=100,verbose=F){
  require("CompQuadForm")
  e=AD.weibull.eigen(neig)
  plb=pchisq(a/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(a,lambda=e,epsabs = 1e-9,limit=2^7)#play with eps and limit
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
