#' EDF Goodness-of-Fit tests for Exponential Distribution by Bootstrap
#'
#' This function takes in an i.i.d. random sample, use MLE to estimate Exponential
#' parameters, compute probability integral transforms, and computes Cramer-von Mises
#' and Anderson-Darling statistics. P-values are calculated by M bootstrap.
#'
#' @param x random sample
#' @param M number of bootstrap
#'
#' @return gof.exp.bootstrap computes Cramer-von Mises and Anderson-Darling statistics and their P-values by bootstrap.
#' @export
#'
#' @examples
#' x=rexp(1000,1/2)
#' gof.exp.bootstrap(x,M=1000)
#' gof.exp(x)
gof.exp.bootstrap<-function(x, M=10000){
  a2 <- AD.exp(x)
  w2 <- CvM.exp(x)
  n <- length(x)
  pars <- estimate.exp(x)
  scale <- pars
  dat <- rexp(n*M,rate=1/scale)
  dat <- matrix(dat,nrow=M)
  a2vals <- apply(dat,1,AD.exp)
  w2vals <- apply(dat,1,CvM.exp)
  a.pv <- length(a2vals[a2vals>a2])/M
  w.pv <- length(w2vals[w2vals>w2])/M
  w2text <- paste("Cramer-von Mises statistic is ", as.character(round(w2,7)))
  w2text <- paste(w2text,"with P-value is ", as.character(round(w.pv,7)),"\n")
  a2text <- paste("Anderson-Darling statistic is ", as.character(round(a2,7)))
  a2text <- paste(a2text,"with P-value is ", as.character(round(a.pv,7)),"\n")
  cat(w2text)
  cat(a2text)
  invisible(list(Asq = a2,Asq.pvalue =a.pv, Wsq=w2, Wsq.pvalue = w2))
}
