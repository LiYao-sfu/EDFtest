#' AD.laplace
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
AD.laplace = function(x){
  theta = estimate.laplace(x)
  z = cdf.laplace((x-theta[1])/theta[2])
  AD(z)
}
