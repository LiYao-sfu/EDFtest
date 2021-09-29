#' P-value of EDF statistics W^2 for Gamma Distribution
#'
#' Compute p-value of the given Cramer-von Mises statistic W^2
#'
#' @param w W^2 for Gamma Distribution
#' @param shape shape parameter of Gamma distribution
#' @param neig number of eigenvalues
#' @param verbose logical; if TRUE, print warning messages
#'
#' @return CvM.gamma.pvalue gives p-value of the Cramer-von Mises statistic of a uniform sample.
#' @export
#'
#' @examples
#' x = rgamma(100,1,1)
#' wsq = CvM.gamma(x)
#' CvM.gamma.pvalue(wsq,1)
CvM.gamma.pvalue = function(w,shape , neig = 100,verbose=FALSE){
  require("CompQuadForm")
  e = CvM.gamma.eigen(neig,shape=shape)
  plb=pchisq(w/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(w,lambda=e)#play with eps and limit
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

