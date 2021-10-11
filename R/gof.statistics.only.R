#' Goodness-of-Fit tests for statistics only
#'
#' @description
#' This function takes in a vector of probability integral transforms
#' and computes whichever of the three statistics is asked for.
#' By default {A^2},{W^2}, and {U^2} are computed
#'
#' @param pit A vector of probability integral transforms.
#' @param AD Logical; if TRUE print Anderson-Darling statistic.
#' @param CvM Logical; if TRUE print Cramér-von Mises statistic.
#' @param Watson Logical; if TRUE print Watson statistic.
#'
#' @return `gof.statistics.only` gives Anderson-Darling, Cramér-von Mises or Watson statistics
#' of a given probability integral transforms
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
