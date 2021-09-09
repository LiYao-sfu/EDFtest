#' EDF Goodness-of-Fit tests for Statistics Only
#'
#' This function takes in a vector of probability integral transforms
#' and computes whichever of the three statistics is asked for.
#' By default A^2,W^2, and U^2 are computed
#'
#' @param pit a vector of probability integral transforms
#' @param a2 logical, if TRUE print Anderson-Darling statistic
#' @param w2 logical, if TRUE print Cramer-von Mises statistic
#' @param u2 logical, if TRUE print Watson statistic
#'
#' @return gof.statistics.only gives Anderson-Darling, Cramer-von Mises or Watson statistics
#' @export
#'
#' @examples
#' x= runif(100)
#' gof.statistics.only(x)
gof.statistics.only=function(pit,a2=TRUE,w2=TRUE,u2=TRUE){
  if(any(pit<0))stop("Negative Probability Integral Transforms Not Allowed")
  if(any(pit>1))stop("Probability Integral Transforms More than 1 Not Allowed")
  p=sort(pit)
  out=list()
  if(a2){out$A2 = AD(p)}
  if(w2){out$W2 = CvM(p)}
  if(u2){out$U2 = Watson(p)}
  return(out)
}
