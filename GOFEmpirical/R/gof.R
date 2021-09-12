#' EDF Goodness-of-Fit tests
#'
#' This function calculate goodness-of-fit test statistics and their P-values.
#' The two statistics computed are the Empirical Distribution function statistics
#' called Cramér-von Mises and Anderson-Darling statistic.
#'
#' The statistics and their P-values can be used to assess an assumed distribution.
#' You have an iid sample from some distribution $F$ and want to test the hypothesis
#' that $F$ is a member of some specific parametric family. The following families are available:
#'
#' Normal(\mu,\sigma^2) Gamma(shape = \alpha, scale = \beta)
#' Logistic(location = \mu, scale = \beta) Laplace(location = \mu, scale = \beta)
#' Weibull(shape = \alpha, scale = \beta) Extreme Value(location = \mu, scale = \beta)
#'
#' @param x random sample
#' @param family assumed distribution
#' @param print Logical, if TRUE print statistics and P-values
#' @param verbose Logical, if TRUE print every step
#' @param bootstrap Logical, if TRUE use bootstrap method, otherwise use Imhof
#' @param M number of bootstrap
#'
#' @return gof computes Cramer-von Mises and Anderson-Darling statistics and their P-values.
#' @export
#'
#' @examples
#' x1=rexp(1000)
#' gof(x=x1,family="Exponential",print=TRUE,verbose=FALSE,bootstrap=FALSE)
#' gof(x=x1,family="Exponential",print=TRUE,verbose=FALSE,bootstrap=TRUE,M=1000)
#'
#' x2=rnorm(1000)
#' gof(x=x2,family="Normal",print=TRUE,verbose=FALSE,bootstrap=FALSE)
#' gof(x=x2,family="Normal",print=TRUE,verbose=FALSE,bootstrap=TRUE,M=1000)
#'
#' x3=rgamma(1000,2)
#' gof(x=x3,family="Gamma",print=TRUE,verbose=FALSE,bootstrap=FALSE)
#' gof(x=x3,family="Gamma",print=TRUE,verbose=FALSE,bootstrap=TRUE,M=1000)
#'
#' x4=rweibull(1000,2)
#' gof(x=x4,family="Weibull",print=TRUE,verbose=FALSE,bootstrap=FALSE)
#' gof(x=x4,family="Weibull",print=TRUE,verbose=FALSE,bootstrap=TRUE,M=1000)
#'
#' x5=rlogis(1000)
#' gof(x=x5,family="Logistic",print=TRUE,verbose=FALSE,bootstrap=FALSE)
#' gof(x=x5,family="Logistic",print=TRUE,verbose=FALSE,bootstrap=TRUE,M=1000)
#'
#' library(L1pack)
#' x6=rlaplace(1000)
#' gof(x=x6,family="Laplace",print=TRUE,verbose=FALSE,bootstrap=FALSE)
#' gof(x=x6,family="Laplace",print=TRUE,verbose=FALSE,bootstrap=TRUE,M=1000)
gof=function(x,family,print=TRUE,verbose=FALSE,bootstrap=FALSE,M=10000){
  if(!bootstrap){method="Imhof"} else{method='Bootstrap'}

  cat("Method:",method)
  cat(", assumed distribution:",family,"\n")

  if(family=="Exponential"& !bootstrap){gof.exp(x,print=print,verbose=verbose)}
  else if(family=="Exponential"& bootstrap){gof.exp.bootstrap(x,M=M)}

  else if(family=="Gamma"& !bootstrap){gof.gamma(x,print=print,verbose=verbose)}
  else if(family=="Gamma"& bootstrap){gof.gamma.bootstrap(x,M=M)}

  else if(family=="Laplace"& !bootstrap){gof.laplace(x,print=print,verbose=verbose)}
  else if(family=="Laplace"& bootstrap){gof.laplace.bootstrap(x,M=M)}

  else if(family=="Logistic"& !bootstrap){gof.logistic(x,print=print,verbose=verbose)}
  else if(family=="Logistic"& bootstrap){gof.logistic.bootstrap(x,M=M)}

  else if(family=="Normal"& !bootstrap){gof.normal(x,print=print,verbose=verbose)}
  else if(family=="Normal"& bootstrap){gof.normal.bootstrap(x,M=M)}

  else if(family=="Weibull"& !bootstrap){gof.weibull(x,print=print,verbose=verbose)}
  else if(family=="Weibull"& bootstrap){gof.weibull.bootstrap(x,M=M)}

  else{warning('Please enter valid parameters')}
}
