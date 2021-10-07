#' EDF Statistic W^2 for a Given Distribution
#'
#' @description
#' Compute Cramer-von Mises goodness-of-fit statistic W^2 for an iid sample, x, to test for
#' the given distribution with parameters unknown. Estimate parameters by ML
#' using \code{EDFtest} mle function by default.
#'
#' @details
#'
#' @param x A random sample.
#' @param parameter Parameters of a given distribution.
#'
#' @return Cramer-von Mises statistic of a uniform sample.
#' @export
#'
#' @examples
#' x0=runif(n=100,min=-1,max=1)
#' CvM.uniform(x0)
#'
#' x1=rnorm(n=100,mean=0,sd=1)
#' CvM.normal(x1)
#'
#' x2=rgamma(n=100,shape=1,scale=1)
#' CvM.gamma(x2)
#'
#' x3=rlogis(n=100,location=0,scale=1)
#' CvM.logistic(x3)
#'
#' x4= rmutil::rlaplace(n=100,m=0,s=1)
#' CvM.laplace(x4)
#'
#' x5=rweibull(n=100,shape=1,scale=1)
#' CvM.weibull(x5)
#'
#' x6=rexp(n=100,rate=1/2)
#' CvM.exp(x6)
CvM.uniform = function(x,parameter=estimate.uniform(x)){
  z <- cdf.uniform(x,parameter)
  CvM(z)
}

#' @export
#' @rdname CvM.uniform
CvM.normal = function(x,parameter=estimate.normal(x)){
  z <- cdf.normal(x,parameter)
  CvM(z)
}

#' @export
#' @rdname CvM.uniform
CvM.gamma <- function(x,parameter=estimate.gamma(x)){
  z <- cdf.gamma(x,parameter)
  CvM(z)
}

#' @export
#' @rdname CvM.uniform
CvM.logistic <- function(x,parameter=estimate.logistic(x)){
  z <- cdf.logistic(x,parameter)
  CvM(z)
}

#' @export
#' @rdname CvM.uniform
CvM.laplace = function(x,parameter=estimate.laplace(x)){
  z = cdf.laplace(x,parameter)
  CvM(z)
}

#' @export
#' @rdname CvM.uniform
CvM.weibull <- function(x,parameter=estimate.weibull(x)){
  z <- cdf.weibull(x,parameter)
  CvM(z)
}

#' @export
#' @rdname CvM.uniform
CvM.exp = function(x,parameter=estimate.exp(x)){
  z = cdf.exp(x,parameter)
  CvM(z)
}

# Helpers -----------------------------------------------------------------

CvM <- function(z){
    n <- length(z)
    u <- sort(z)
    i <- (2*(1:n)-1)/(2*n)
    sum((u-i)^2)+1/(12*n)
}

# Compute the probability integral transforms of
# the random sample x in a given distribution
cdf.uniform = function(x,theta){
  punif(x,min=theta[1],max=theta[2])
}

cdf.normal = function(x,theta){
  pnorm(x,mean=theta[1],sd=theta[2])
}

cdf.gamma = function(x,theta){
  pgamma(x,shape=theta[1],scale=theta[2])
}

cdf.logistic = function(x,theta){
  plogis(x,location=theta[1],scale=theta[2])
}

cdf.laplace = function(x,theta){
  m = theta[1]
  s = theta[2]
  v = exp(x-m)/2
  v[x>m] = 1-exp(-(x-m)[x>m])/2
  v
}

cdf.weibull = function(x,theta){
  pweibull(x,shape=theta[1],scale=theta[2])
}

cdf.exp = function(x,theta){
  pexp(x,rate=1/theta)
}







