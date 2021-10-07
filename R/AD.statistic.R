#' EDF Statistic A^2 for a Given Distribution
#'
#' @description
#' Compute the Anderson-Darling goodness-of-fit statistic A^2 for an iid sample,
#' x, to test for the given distribution with parameters unknown.
#' Estimate parameters by ML using \code{EDFtest} mle function by default.
#'
#' @details
#'
#' @param x A random sample.
#' @param parameter Parameters of a given distribution.
#'
#' @return Anderson-Darling statistic of a uniform sample.
#' @export
#'
#' @examples
#' x0=runif(n=100,min=-1,max=1)
#' AD.uniform(x0)
#'
#' x1=rnorm(n=100,mean=0,sd=1)
#' AD.normal(x1)
#'
#' x2=rgamma(n=100,shape=1,scale=1)
#' AD.gamma(x2)
#'
#' x3=rlogis(n=100,location=0,scale=1)
#' AD.logistic(x3)
#'
#' x4= rmutil::rlaplace(n=100,m=0,s=1)
#' AD.laplace(x4)
#'
#' x5=rweibull(n=100,shape=1,scale=1)
#' AD.weibull(x5)
#'
#' x6=rexp(n=100,rate=1/2)
#' AD.exp(x6)
AD.uniform = function(x,parameter=estimate.uniform(x)){
  s <- sort(x)[-c(1,length(x))]
  z <- cdf.uniform(s,parameter)
  return(AD(z))
}

#' @export
#' @rdname AD.uniform
AD.normal = function(x,parameter=estimate.normal(x)){
  z <- cdf.normal(x,parameter)
  AD(z)
}

#' @export
#' @rdname AD.uniform
AD.gamma <- function(x,parameter=estimate.gamma(x)){
  z <- cdf.gamma(x,parameter)
  AD(z)
}

#' @export
#' @rdname AD.uniform
AD.logistic <- function(x,parameter=estimate.logistic(x)){
  z <- cdf.logistic(x,parameter)
  AD(z)
}

#' @export
#' @rdname AD.uniform
AD.laplace = function(x,parameter=estimate.laplace(x)){
  z = cdf.laplace(x,parameter)
  AD(z)
}

#' @export
#' @rdname AD.uniform
AD.weibull <- function(x,parameter=estimate.weibull(x)){
  z <- cdf.weibull(x,parameter)
  AD(z)
}

#' @export
#' @rdname AD.uniform
AD.exp = function(x,parameter=estimate.exp(x)){
  z = cdf.exp(x,parameter)
  AD(z)
}

AD <- function(z){
  #AKA AD.uniform()
  n <- length(z)
  u <- sort(z)
  i <- 2*(1:n)-1
  -n -sum(i*(log(u)+log(1-rev(u))))/n
}
