#' P-value of EDF statistics A^2 for Gamma Distribution
#'
#' Compute p-value of the given Anderson-Darling statistic A^2
#'
#' @param a A^2 for Gamma Distribution
#' @param shape shape parameter for Gamma distribution
#' @param neig number of eigenvalues
#' @param verbose logical; if TRUE, print warning messages
#'
#' @return AD.gamma.pvalue gives p-value of the Anderson-Darling statistic of a uniform sample.
#' @export
#'
#' @examples
#' x = rgamma(100,1,2)
#' asq = AD.gamma(x)
#' AD.gamma.pvalue(asq,1)
AD.gamma.pvalue = function(a,shape,neig = 100,verbose=FALSE){
  require("CompQuadForm")
  e = AD.gamma.eigen(neig,shape=shape)
  plb=pchisq(a/max(e),df=1,lower.tail = FALSE)
  warn=getOption("warn")
  im = imhof(a,lambda=e)#play with eps and limit
  options(warn=warn)
  aerror=im$abserr
  p=im$Qq
  if(p<0&&verbose)cat("for A = ",a," and neig = ",neig,
                      " imhof returned a negative probability\n")
  if(p<plb){
    p=plb
    if(verbose) cat("forA = ",a," and neig = ",neig,
                    " p was replaced by a lower bound on p: ",plb, "\n")
  }
  list(P=p,error=aerror)
}
