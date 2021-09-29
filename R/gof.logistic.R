#' EDF Goodness-of-Fit tests for Logistic Distribution
#'
#' This function takes in an i.i.d. random sample, use MLE to estimate Logistic
#' parameters, compute probability integral transforms, and computes Cramer-von Mises
#' and Anderson-Darling statistics and their P-values
#'
#' @param x random sample
#' @param print Logical, if TRUE print statistics and P-values
#' @param verbose Logical, if TRUE print every step
#'
#' @return gof.logistic computes Cramer-von Mises and Anderson-Darling statistics and their P-values.
#' @export
#'
#' @examples
#' x = rlogis(100,1)
#' gof.logistic(x)
#' gof.logistic(x,print=TRUE,verbose=TRUE)
gof.logistic=function(x,print=T,verbose=F){
  require("CompQuadForm")

  #  Estimate the parameters
  pars=estimate.logistic(x)
  if(verbose){cat("log-Logistic parameter estimates", pars, "\n")}

  #  Compute the pit
  pit=plogis(x,location=pars[1],scale=pars[2])
  if(verbose){cat("PITs are done \n \n")}

  #  Compute two gof statistics
  w = CvM(pit)
  a = AD(pit)
  if(verbose){cat("Statistics are ",w,a," \n \n")}

  #  Compute their p-values
  w.p=CvM.logistic.pvalue(w,verbose=verbose)
  if(verbose){
    cat("Cramer von Mises P value output \n")
    print(w.p)
    cat("\n\n")
  }
  w.p=w.p$P
  a.p=AD.logistic.pvalue(a)$P
  if(verbose){
    cat("Anderson-Darling P value output \n")
    print(a.p)
    cat("\n\n")
  }
  if(print){
    cat("Cramer-von Mises statistic is ",w,"with P-value is ",w.p,"\n")
    cat("Anderson-Darling statistic is ",a,"with P-value is ",a.p,"\n")
  }
  invisible(list(w=w,w.p=w.p,a=a,a.p=a.p))
}
