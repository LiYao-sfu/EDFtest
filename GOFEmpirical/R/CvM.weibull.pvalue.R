#' CvM.weibull.pvalue
#'
#' @param w
#' @param neig
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
CvM.weibull.pvalue = function(w,neig=100,verbose=FALSE){
  library("CompQuadForm")
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
