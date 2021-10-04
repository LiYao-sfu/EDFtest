#' EDF statistics W^2 for Poisson Distribution
#'
#' Compute Cramer-von Mises statistic W^2 for an iid sample, x, to test for the Poisson distribution with parameters unknown.
#' Estimate parameters by ML using "estimate.poisson" by default.
#'
#' @param x random sample
#' @param eps how much prob you are willing to omit
#'
#' @return CvM.poisson gives Cramer-von Mises statistic of a uniform sample.
#' @export
#'
#' @examples
#' x= rpois(100,2)
#' CvM.poisson(x)
CvM.poisson = function(x,eps=10^(-9)){
  n=length(x)
  s=mean(x)
  obs=tabulate(x+1)
  #    print(obs)
  m=max(x)
  p=dpois(0:m,s)
  tw=p
  #
  # The Cramer von Mises statistic uses weights corresponding to dF
  #  in the continuous case.  There are at least 2 reasonable possible
  #  formulas for the weights.  We use our original choice here
  #
  # The formulas below are for weights which would give the same statistic if the #  cells were reversed in order.
  #
  #    p1=dpois(1:(m+1),s)
  #
  #    tw=(p+p1)/2
  #
  e = n*p
  z=cumsum(obs-e)
  #
  # The statistic is technically an infinite sum which we truncate
  #  in the code below to miss only eps of the tail probability.
  #  If the more elaborate weights above are used the code below
  #  needs to change, too.
  #
  m2=max(qpois(1-eps,s),100)
  p2=dpois((m+1):m2,s)
  q2=ppois(m:(m2-1),s,lower.tail=F)
  tail=n*sum(q2^2*p2)
  sum(z^2*tw)/n+tail
}

#' P-value of EDF statistics W^2 for Poisson Distribution by bootstrap
#'
#' Compute the Cramer-von Mises statistic W^2 and a corresponding P-value
#' for testing the hypothesis that a sample comes from a Poisson distribution
#' with mean to be estimated from the data
#
#' The P-value is computed by conditional sampling using nmc samples from
#' the conditional distribution of x given sum(x) for a Poisson distribution.
#' This distribution is multinomial.
#'
#' @param x random sample
#' @param nmc
#' @param return.samples
#' @param print logical, if True print the statistic and p-value
#'
#' @return CvM.poisson.pvalue.bootstrap gives Cramer-von Mises statistic and its p-value.
#' @export
#'
#' @examples
#' x= rpois(100,2)
#' CvM.poisson(x)
#' CvM.poisson.pvalue.bootstrap(x,return.samples=T)
CvM.poisson.pvalue.bootstrap = function(x,nmc=200,return.samples=FALSE,print=TRUE){
  w=CvM.poisson(x)
  v= rep(0,nmc)
  n=length(x)
  t.suff=sum(x)
  #    print(c(w,v,n,t.suff))
  u=rmultinom(nmc,t.suff,rep(1,n)/n)
  v=apply(u,2,CvM.poisson)
  p.value=length(v[v >w])/nmc
  if(return.samples) return(list(cvm=w,p.value=p.value,samples=v))
  if(print){
    cat(paste("Cramer-von Mises Statistic ",as.character(round(w,4)),"  "))
    cat(paste("Corresponding P-Value ",as.character(p.value)))
  }
  invisible(list(cvm=w,p.value=p.value))
}
