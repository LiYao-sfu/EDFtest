#' EDF Goodness-of-Fit tests for Gamma Distribution
#'
#' This function takes in an i.i.d. random sample, use MLE to estimate Gamma
#' parameters, compute probability integral transforms, and computes Cramer-von Mises
#' and Anderson-Darling statistics and their P-values
#'
#' @param x random sample
#' @param print Logical, if TRUE print statistics and P-values
#' @param verbose Logical, if TRUE print every step
#'
#' @return gof.gamma computes Cramer-von Mises and Anderson-Darling statistics and their P-values.
#' @export
#'
#' @examples
#' x = rgamma(100,1)
#' gof.gamma(x)
#' gof.gamma(x,print=TRUE,verbose=TRUE)
gof.gamma=function(x,print=TRUE,verbose=FALSE){
  require("CompQuadForm")

  #  Estimate the parameters
  pars=estimate.gamma(x)
  if(verbose){cat("Gamma parameter estimates", pars, "\n")}

  #  Compute the pit
  pit=pgamma(x,shape=pars[1],scale=pars[2])
  if(verbose){cat("PITs are done \n \n")}

  #  Compute two gof statistics
  w = CvM(pit)
  a = AD(pit)

  #  Compute their p-values
  w.p=CvM.gamma.pvalue(w,shape=pars[1])$P
  if(verbose){
    cat("Cramer von Mises P value output \n")
    print(w.p)
    cat("\n\n")
  }
  a.p=AD.gamma.pvalue(a,shape=pars[1])$P
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
