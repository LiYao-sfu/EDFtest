#' EDF Goodness-of-Fit tests for Exponential Distribution
#'
#' This function takes in an i.i.d. random sample, use MLE to estimate Exponential
#' parameters, compute probability integral transforms, and computes Cramer-von Mises
#' and Anderson-Darling statistics and their P-values
#'
#' @param x random sample
#' @param print Logical, if TRUE print statistics and P-values
#' @param verbose Logical, if TRUE print every step
#'
#' @return gof.exp computes Cramer-von Mises and Anderson-Darling statistics and their P-values.
#' @export
#'
#' @examples
#' x = rexp(100,1/2)
#' gof.exp(x)
#' gof.exp(x,print=TRUE,verbose=TRUE)
gof.exp=function(x,print=TRUE,verbose=FALSE){
  require("CompQuadForm")

  #  Estimate the parameters
  pars=estimate.exp(x)
  if(verbose){cat("Exponential parameter estimates", pars, "\n")}

  #  Compute the pit
  pit=pexp(x,rate=1/pars)
  if(verbose){cat("PITs are done \n \n")}

  #  Compute two gof statistics
  w = CvM(pit)
  a = AD(pit)

  #  Compute their p-values
  w.p=CvM.exp.pvalue(w)$P
  if(verbose){
    cat("Cramer von Mises P value output \n")
    print(w.p)
    cat("\n\n")
  }
  a.p=AD.exp.pvalue(a)$P
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
