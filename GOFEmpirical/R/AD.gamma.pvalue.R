#' AD.gamma.pvalue
#'
#' @param a
#' @param shape
#' @param neig
#' @param verbose
#'
AD.gamma.pvalue = function(a,shape,neig = 100,verbose=FALSE){
  library("CompQuadForm")
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
