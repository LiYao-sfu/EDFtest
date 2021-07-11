#' AD.exp
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
AD.exp = function(x){
  theta = estimate.exp(x)
  z = pexp(x/theta[1],rate=1)
  AD(z)
}
