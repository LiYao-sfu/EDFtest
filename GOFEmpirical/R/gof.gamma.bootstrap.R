#' EDF Goodness-of-Fit tests for Gamma Distribution by Bootstrap
#'
#' This function takes in an i.i.d. random sample, use MLE to estimate Gamma
#' parameters, compute probability integral transforms, and computes Cramer-von Mises
#' and Anderson-Darling statistics. P-values are calculated by M bootstrap.
#'
#' @param x random sample
#' @param M number of bootstrap
#'
#' @return gof.gamma.bootstrap computes Cramer-von Mises and Anderson-Darling statistics and their P-values by bootstrap.
#' @export
#'
#' @examples
#' x=rgamma(10,1)
#' gof.gamma.bootstrap(x,M=1000)
#' gof.gamma(x)
gof.gamma.bootstrap<-function(x, M=10000){
  a2 <- AD.gamma(x)
  w2 <- CvM.gamma(x)
  n <- length(x)
  pars <- estimate.gamma(x)
  alpha <- pars[1]
  beta <- pars[2]
  dat <- rgamma(n*M,shape=alpha,scale=beta)
  dat <- matrix(dat,nrow=M)

  # Make it efficient!

  # ests <- apply(dat,1,estimate.gamma)
  # pit <- lapply(dat,1,cdf.gamma,theta=ests)
  # dim(pit)
  # cdf.gamma(dat[1,],theta=ests[,1])
  # cdf.gamma(dat[2,],theta=ests[,2])
  # cdf.gamma(dat[2,],theta=ests[,1])
  # pit[,1:2]

  a2vals <- apply(dat,1,AD.gamma)
  w2vals <- apply(dat,1,CvM.gamma)
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
