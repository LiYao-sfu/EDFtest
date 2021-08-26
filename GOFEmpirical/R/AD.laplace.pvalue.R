#' P-value of EDF statistics A^2 for Laplace Distribution
#'
#' Compute p-value of the given Anderson-Darling statistic A^2
#'
#' @param a A^2 for Exponential Distribution
#' @param neig number of eigenvalues
#' @param verbose logical; if TRUE, print warning messages
#'
#' @return AD.laplace.pvalue gives p-value of the Anderson-Darling statistic of a uniform sample.
#' @export
#'
#' @examples
#' library(L1pack)
#' x= rlaplace(100)
#' asq = AD.laplace(x)
#' AD.laplace.pvalue(asq)
AD.laplace.pvalue = function(a,neig=100,verbose=FALSE){
  require("CompQuadForm")
  e = AD.laplace.eigen(neig)
  plb=pchisq(a/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(a,lambda=e,epsabs = 1e-9,limit=2^15)
  options(warn=warn)
  aerror=im$abserr
  p=im$Qq
  if(p<0&&verbose)cat("for A = ",a," and neig = ",neig,
                      " imhof returned a negative probability\n")
  if(p<plb){
    p=plb
    if(verbose) cat("for A = ",a," and neig = ",neig,
                    " p was replaced by a lower bound on p: ",plb, "\n")
  }
  list(P=p,error=aerror)
}
