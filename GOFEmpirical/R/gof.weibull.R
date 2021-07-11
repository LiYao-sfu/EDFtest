#' gof.weibull
#'
#' @param x
#' @param print
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
gof.weibull=function(x,print=T,verbose=F){
  library("CompQuadForm")
  #
  #  Does the weibull distribution
  #
  #  I don't check to see if x is valid input
  #
  #  Estimate the parameters
  pars=estimate.weibull(x)
  if(verbose){
    cat("Weibull parameter estimates", pars, "\n")
  }

  #
  #  Compute the pit
  #
  pit=pweibull(x,shape=pars[1],scale=pars[2])
  if(verbose){
    cat("PITs are done \n \n")
  }
  #
  #
  #  Compute two gof statistics
  #
  w = CvM(pit)
  a = AD(pit)
  #
  #  Now use the large sample theory and the method of imhof
  #
  w.p=CvM.weibull.pvalue(w)$P
  if(verbose){
    cat("Cramer von Mises P value output \n")
    print(w.p)
    cat("\n\n")
  }
  a.p=AD.weibull.pvalue(a)$P
  if(print){
    cat("Cramer-von Mises statistic is ",w," with P = ",w.p,"\n")
    cat("Anderson-Darling statistic is ",a," with P = ",a.p,"\n")
  }
  invisible(list(w=w,w.p=w.p,a=a,a.p=a.p))
}
