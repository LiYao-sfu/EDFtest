#' CvM.exp
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
CvM.exp = function(x){
  z = pexp(x/mean(x),rate=1)
  CvM(z)
}
