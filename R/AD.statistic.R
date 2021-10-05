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
#' x1=rnorm(n=100,mean=0,sd=1)
#' AD.normal(x1)
#' AD.normal(x1,c(0,1))
#'
#' x2=rgamma(n=100,shape=1,scale=1)
#' AD.gamma(x2)
#' AD.gamma(x2,parameter=c(1,1))
#'
#' x3=rlogis(n=100,location=0,scale=1)
#' AD.logistic(x3)
#' AD.logistic(x3,c(0,1))
#'
#' x4= L1pack::rlaplace(n=100,location=0,scale=1)
#' AD.laplace(x4)
#' AD.laplace(x4,c(0,1))
#'
#' x5=rweibull(n=100,shape=1,scale=1)
#' AD.weibull(x5)
#' AD.weibull(x5,c(1,1))
#'
#' x6=rexp(n=100,rate=1/2)
#' AD.exp(x6)
#' AD.exp(x6,parameter=2)
AD.uniform <- function(x,parameter=estimate.uniform(x)){
  z <- cdf.uniform(x,parameter)
  AD(z)
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
