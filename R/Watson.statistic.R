#' EDF Statistic U^2 for a Given Distribution
#'
#' @description
#' Compute the Watson goodness-of-fit statistic U^2 for an iid sample,
#' x, to test for the given distribution with parameters unknown.
#' Estimate parameters by ML using \code{EDFtest} mle function by default.
#'
#' @details
#'
#' @param x A random sample.
#' @param parameter Parameters of a given distribution.
#'
#' @return Watson's statistic of a uniform sample.
#' @export
#'
#' @examples
#' x1=rnorm(n=100,mean=0,sd=1)
#' Watson.normal(x1)
#' Watson.normal(x1,parameter=c(0,1))
#'
#' x2=rgamma(n=100,shape=1,scale=1)
#' Watson.gamma(x2)
#' Watson.gamma(x2,parameter=c(1,1))
#'
#' x3=rlogis(n=100,location=0,scale=1)
#' Watson.logistic(x3)
#' Watson.logistic(x3,parameter=c(0,1))
#'
#' x4= rmutil::rlaplace(n=100,m=0,s=1)
#' Watson.laplace(x4)
#' Watson.laplace(x4,parameter=c(0,1))
#'
#' x5=rweibull(n=100,shape=1,scale=1)
#' Watson.weibull(x5)
#' Watson.weibull(x5,parameter=c(1,1))
#'
#' x6=rexp(n=100,rate=1/2)
#' Watson.exp(x6)
#' Watson.exp(x6,parameter=2)
Watson.normal = function(x,parameter=estimate.normal(x)){
    z <- cdf.normal(x, theta=parameter)
    Watson(z)
}

#' @export
#' @rdname Watson.normal
Watson.gamma = function(x,parameter=estimate.gamma(x)){
    z <- cdf.gamma(x,parameter)
    Watson(z)
}

#' @export
#' @rdname Watson.normal
Watson.logistic = function(x,parameter=estimate.logistic(x)){
    z <- cdf.logistic(x,parameter)
    Watson(z)
}

#' @export
#' @rdname Watson.normal
Watson.laplace = function(x,parameter=estimate.laplace(x)){
    z <- cdf.laplace(x,parameter)
    Watson(z)
}

#' @export
#' @rdname Watson.normal
Watson.weibull = function(x,parameter=estimate.weibull(x)){
    z <- cdf.weibull(x,parameter)
    Watson(z)
}

#' @export
#' @rdname Watson.normal
Watson.exp = function(x,parameter=estimate.exp(x)){
    z <- cdf.exp(x,parameter)
    Watson(z)
}

Watson <- function(z){
    # Given probability integral transforms, compute the Watson goodness-of-fit statistic
    # for testing uniformity on the unit interval.
    n <- length(z)
    u <- sort(z)
    i <- (2*(1:n)-1)/(2*n)
    w=sum((u-i)^2)+1/(12*n)
    corr=n*(mean(u)-0.5)^2 #correction
    w-corr
}
