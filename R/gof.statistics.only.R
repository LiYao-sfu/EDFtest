#' EDF Goodness-of-Fit tests for Statistics Only
#'
#' @description
#' This function takes in a vector of probability integral transforms
#' and computes whichever of the three statistics is asked for.
#' By default A^2,W^2, and U^2 are computed
#'
#' @details
#'
#' @param pit A vector of probability integral transforms.
#' @param a2 Logical; if TRUE print Anderson-Darling statistic
#' @param w2 Logical; if TRUE print Cramer-von Mises statistic
#' @param u2 Logical; if TRUE print Watson statistic
#'
#' @return `gof.statistics.only` gives Anderson-Darling, Cramer-von Mises or Watson statistics
#' @export
#'
#' @examples
#' x= runif(100)
#' gof.statistics.only(x)
gof.statistics.only=function(pit,AD=TRUE,CvM=TRUE,Watson=TRUE){
  if(any(pit<0))stop("Negative Probability Integral Transforms Not Allowed")
  if(any(pit>1))stop("Probability Integral Transforms More than 1 Not Allowed")
  p=sort(pit)
  out=list()
  if(AD){out$AD = AD(p)}
  if(CvM){out$CvM = CvM(p)}
  if(Watson){out$Watson = Watson(p)}
  return(out)
}
