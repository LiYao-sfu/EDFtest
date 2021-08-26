#' MLE for Weibull Distribution
#'
#' Estimate shape and scale parameters of the Weibull distribution by the method of maximum likelihood.
#'
#' @param x random sample
#'
#' @return estimated shape and scale parameters of the Weibull distribution.
#' @export
#'
#' @examples
#' x=rweibull(10,1)
#' estimate.weibull(x)
estimate.weibull <- function(x){
    n <- length(x)
    m1 <- mean(x)
    m2 <- var(x)
    b <- m2/m1
    a <- m1/b
    mlog <- mean(log(x))
    aold <- a
    ml <- mean(log(x))
    mla <- mean(x^aold*log(x))
    ma <- mean(x^aold)
    anew <- 1/(mla/ma-ml)
    anew <- sqrt(anew*aold)

    while ( abs(anew-aold) > 1e-7){
      aold <- anew
      mla <- mean(x^aold*log(x))
      ma <- mean(x^aold)
      anew <- 1/(mla/ma-ml)
      anew <- sqrt(anew*aold)
    }
    beta <- mean(x^anew)^(1/anew)
    alpha <- anew
    c(alpha, beta)
  }
