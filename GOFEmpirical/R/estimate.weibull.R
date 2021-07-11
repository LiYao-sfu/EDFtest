#' MLE for Weibull Distribution
#'
#'
#'
#' @param x numerical vector
#'
#' @return vector of two parameters, the first alpha, the second beta
#' @export
#'
#' @examples
#'
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
