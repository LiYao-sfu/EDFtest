#' AD.normal.eigen
#'
#' @param n
#'
#' @return
#' @export
#'
#' @examples
AD.normal.eigen = function(n){
  # mean.asq.weibull = .4519444511 # from Maple
  M=CvM.normal.covmat(n)
  e=eigen(M)$values/n
  e # *mean.asq.weibull/sum(e)
}
