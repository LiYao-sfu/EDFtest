#' CvM.weibull.pvalue
#'
#' @param w
#' @param nd
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
CvM.weibull.pvalue = function(w,nd=9,verbose=F){
  library("CompQuadForm")
  p=rep(0,nd)
  aerror=p
  ierror=0
  emat = matrix(0,ncol=2^nd,nrow=nd)
  for(i in nd:nd){
    emat[i,1:2^i]=CvM.weibull.eigen(2^i)
    plb=pchisq(w/max(emat[i,]),df=1,lower.tail = FALSE)
    warn=getOption("warn")
    im = imhof(w,lambda=emat[i,],epsabs = 1e-9,limit=2^15)
    options(warn=warn)
    aerror[i]=im$abserr
    p[i]=im$Qq
    if(p[i]<aerror[i]){
      ierror=1
      if(p[i]<0&&verbose)cat("for i = ",i," of nd = ",nd," imhof returned a negative probability\n")
    }
    if(p[i]<plb){
      p[i]=plb
      if(verbose)cat("for i = ",i," of nd = ",nd," p =",p[i]," was replaced by a lower bound on p:",plb, "\n")
    }
  }
  if(ierror >0&&verbose){
    cat("Large errors detected\n")
    cat("Here is the table of computed p-values before extrapolation\n")
    cat(p)
    cat("\n")
  }
  #p.ext.1=2*p[nd]-p[nd-1]
  p.ext.1=p[nd]
  if(abs(p.ext.1-p[nd])>0.00001&&verbose){
    ierror=1
    cat("Richardson extrapolation changed P by more than 0.00001\n")
  }

  list(P=p.ext.1,all.P=p,error=aerror)
}
