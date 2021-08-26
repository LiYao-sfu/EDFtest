#' P-value of EDF statistics W^2 for Laplace Distribution
#'
#' Compute p-value of the given Cramer-von Mises statistic W^2
#'
#' @param w W^2 for Laplace Distribution
#' @param neig number of eigenvalues
#' @param verbose logical; if TRUE, print warning messages
#'
#' @return CvM.laplace.pvalue gives p-value of the Cramer-von Mises statistic of a uniform sample.
#' @export
#'
#' @examples
#' library(L1pack)
#' x = rlaplace(1000,0,1)
#' wsq = CvM.laplace(x)
#' CvM.laplace.pvalue(wsq)
CvM.laplace.pvalue = function(w,neig=100,verbose=F){
  require("CompQuadForm")
  e = CvM.laplace.eigen(neig)
  plb=pchisq(w/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(w,lambda=e,epsabs = 1e-9,limit=2^15)
  options(warn=warn)
  aerror=im$abserr
  p=im$Qq
  if(p<0&&verbose)cat("for w = ",w," and neig = ",neig,
                      " imhof returned a negative probability\n")
  if(p<plb){
    p=plb
    if(verbose) cat("for w = ",w," and neig = ",neig,
                    " p was replaced by a lower bound on p: ",plb, "\n")
  }
  list(P=p,error=aerror)
}
