#' P-value of EDF statistics W^2 for Normal Distribution
#'
#' Compute p-value of the given Cramer-von Mises statistic W^2
#'
#' @param w W^2 for Normal Distribution
#' @param neig number of eigenvalues
#' @param verbose logical; if TRUE, print warning messages
#'
#' @return CvM.normal.pvalue gives p-value of the Cramer-von Mises statistic of a uniform sample.
#' @export
#'
#' @examples
#' x = rnorm(1000)
#' wsq = CvM.normal(x)
#' CvM.normal.pvalue(wsq)
CvM.normal.pvalue = function(w,neig=100,verbose=F){
  require("CompQuadForm")
  e = CvM.normal.eigen(neig)
  plb=pchisq(w/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(w,lambda=e,epsabs = 1e-9,limit=2^15)
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
