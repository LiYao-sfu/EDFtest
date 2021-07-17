#' CvM.gamma.pvalue
#'
#' @param w
#' @param shape
#' @param neig
#' @param verbose
#'
CvM.gamma.pvalue = function(w,shape , neig = 100,verbose=FALSE){
  library("CompQuadForm")
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

