#' MLE for Gamma Distribution
#'
#' Estimate shape and scale parameters of the Gamma distribution by the method of maximum likelihood.
#' Use the digamma and trigamma functions in Base R and do Newton-Raphson on the profile log-likelihood for alpha
#'
#' @param x random sample
#'
#' @return estimated shape and scale parameters of the Gamma distribution.
#' @export
#'
#' @examples
#' x=rgamma(10,2,1)
#' estimate.gamma(x)
estimate.gamma <- function(x){
  n <- length(x)
  m1 <- mean(x)
  m2 <- var(x)
  b <- m2/m1
  a <- m1/b
  mlog <- mean(log(x))
  logm=log(m1)
  aold <- a
  anew <- aold -(log(aold)-logm + mlog -digamma(aold))/(1/aold-trigamma(aold))
  bnew=m1/anew
  if( anew < 0) anew <- aold/2
  while ( abs(anew-aold) > 1e-7){
    aold <- anew
    old.score = (log(aold)-log(m1)+ mlog -digamma(aold))
    old.score.derivative = 1/aold-trigamma(aold)
    anew <- aold - old.score/old.score.derivative
    if( anew < 0) anew <- aold/2
  }
  beta <- m1/anew
  alpha <- anew
  c(alpha, beta)
}
