#' AD.normal.pvalue
#'
#' @param a
#' @param nd
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
#'
AD.normal.pvalue = function(a,nd=9,verbose=F){
  library("CompQuadForm")
  ierror=0
  p=rep(0,nd)
  aerror=p
  emat = matrix(0,ncol=2^nd,nrow=nd)
  for(i in 1:nd){
    emat[i,1:2^i]=AD.normal.eigen(2^i)
    plb=pchisq(a/max(emat[i,]),df=1,lower.tail = FALSE)
    warn=getOption("warn")
    im = imhof(a,lambda=emat[i,],epsabs = 1e-9,limit=2^7)
    options(warn=warn)
    aerror[i]=im$abserr
    p[i]=im$Qq
    if(p[i]<aerror[i]){
      ierror=1
      if(p[i]<0&&verbose)cat("for i = ",i," of nd = ",nd," imhof returned a negative probability\n")
    }
    if(p[i]<plb){
      p[i]=plb
      if(verbose) cat("for i = ",i," of nd = ",nd," p =",p[i]," was replaced by a lower bound on p:",plb, "\n")
    }
  }
  if(ierror >0 && verbose){
    cat("Large errors detected\n")
    cat("Here is the table of computed p-values before extrapolation\n")
    cat(p)
    cat("\n")
  }
  #  for(j in 2:nd)pext[j,1:(nd+1-j)]=2*pext[j-1,2:(nd+2-j)]-pext[j-1,1:(nd+1-j)]
  p.ext.1=2*p[nd]-p[nd-1]
  if(abs(p.ext.1-p[nd])>0.00001&&verbose){
    ierror=1
    cat("Richardson extraploation changed P by more than 0.00001\n")
  }
  list(P=p.ext.1,all.P=p,error=aerror)
}
