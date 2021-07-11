#' AD.logistic.covmat
#'
#' @param n
#'
#' @return
#' @export
#'
#' @examples
AD.logistic.covmat=function(n){
  s=1:n
  s=s/(n+1)
  t=s
  CvM.logistic.covmat(n)/sqrt(outer(s*(1-s),t*(1-t)))
}
