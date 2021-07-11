#' gof.logistic
#'
#' @param x
#' @param print
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
gof.logistic=function(x,print=T,verbose=F){
  #

  #
  library("CompQuadForm")
  #
  #  Does the logistic distribution
  #  For use with times you must take logs first
  #
  #  I don't check to see if x is valid input
  #
  #  Estimate the parameters
  pars=estimate.logistic(x)
  if(verbose){
    cat("log-Logistic parameter estimates", pars, "\n")
  }
  #
  #  Compute the pit
  #
  pit=plogis(x,location=pars[1],scale=pars[2])
  if(verbose){
    cat("PITs are done \n \n")
  }
  #
  #  Compute two gof statistics
  #
  w = CvM(pit)
  a = AD(pit)
  if(verbose){
    cat("Statistics are ",w,a," \n \n")
  }
  #
  #  Now use the large sample theory and the method of imhof
  #
  w.p=CvM.logistic.pvalue(w,verbose=verbose)
  if(verbose){
    cat("Cramer von Mises P value output \n")
    print(w.p)
    cat("\n\n")
  }
  w.p=w.p$P
  a.p=AD.logistic.pvalue(a)$P
  if(print){
    cat("Cramer-von Mises statistic is ",w," with P = ",w.p,"\n")
    cat("Anderson-Darling statistic is ",a," with P = ",a.p,"\n")
  }
  invisible(list(w=w,w.p=w.p,a=a,a.p=a.p)) #make it a vector
}
