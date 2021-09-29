#' EDF Goodness-of-Fit tests for Laplace Distribution
#'
#' This function takes in an i.i.d. random sample, use MLE to estimate Laplace
#' parameters, compute probability integral transforms, and computes Cramer-von Mises
#' and Anderson-Darling statistics and their P-values
#'
#' @param x random sample
#' @param print Logical, if TRUE print statistics and P-values
#' @param verbose Logical, if TRUE print every step
#'
#' @return gof.laplace computes Cramer-von Mises and Anderson-Darling statistics and their P-values.
#' @export
#'
#' @examples
#' library(rmutil)
#' x= rmutil::rlaplace(1000,0,1)
#' gof.laplace(x)
#' gof.laplace(x,print=TRUE,verbose=TRUE)
gof.laplace=function(x,print=TRUE,verbose=FALSE){
  require("CompQuadForm")
  require("rmutil")

  #  Estimate the parameters
  pars=estimate.laplace(x)
  if(verbose){cat("Laplace parameter estimates", pars, "\n")}

  #  Compute the pit
  pit=rmutil::plaplace(x,m=pars[1],s=pars[2])
  if(verbose){cat("PITs are done \n \n")}

  #  Compute two gof statistics
  w = CvM(pit)
  a = AD(pit)

  #  Compute their p-values
  w.p=CvM.laplace.pvalue(w)$P
  if(verbose){
    cat("Cramer von Mises P value output \n")
    print(w.p)
    cat("\n\n")
  }
  a.p=AD.laplace.pvalue(a)$P
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
