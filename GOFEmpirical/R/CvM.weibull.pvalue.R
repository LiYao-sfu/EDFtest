#' P-value of EDF statistics W^2 for Weibull Distribution
#'
#' Compute p-value of the given Cramer-von Mises statistic W^2
#'
#' @param w W^2 for Weibull Distribution
#' @param neig number of eigenvalues
#' @param verbose logical; if TRUE, print warning messages
#'
#' @return CvM.weibull.pvalue gives p-value of the Cramer-von Mises statistic of a uniform sample.
#' @export
#'
#' @examples
#' x = rweibull(1000,1)
#' wsq = CvM.weibull(x)
#' CvM.weibull.pvalue(wsq)
CvM.weibull.pvalue = function(w,neig=100,verbose=FALSE){
  require("CompQuadForm")
  e=CvM.weibull.eigen(neig)
  plb=pchisq(w/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(w,lambda=e,epsabs = 1e-9,limit=2^7)#play with eps and limit
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
