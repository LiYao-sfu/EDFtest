#' CvM.poisson.pvalue.bootstrap
#'
#' @param x random sample
#' @param nmc
#' @param return.samples
#' @param print
#'
#' @return
#' @export
#'
#' @examples
CvM.poisson.pvalue.bootstrap =
  function(x,nmc=200,return.samples=F,print=TRUE){
    #
    #  This function computes the Cramer-von Mises
    #   statistic and a corresponding P-value
    #   for testing the hypothesis that a sample comes
    #   from a Poisson distribution with mean to be
    #   estimated from the data
    #
    #  The P-value is computed by conditional sampling using
    #   nmc samples from the conditional distribution of x
    #   given sum(x) for a Poisson distribution.  This distribution
    #   is multinomial.
    #
    #  The actual statistic is computed by the function wsq.logistic.pvalue
    #
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
