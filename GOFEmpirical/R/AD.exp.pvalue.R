#' AD.exp.pvalue
#'
#' @param a
#' @param neig
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
AD.exp.pvalue = function(a,neig=100,verbose=FALSE){
  library("CompQuadForm")
  e=AD.exp.eigen(neig)
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
