#' CvM.laplace
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
CvM.laplace = function(x){
  theta = estimate.laplace(x)
  z = cdf.laplace((x-theta[1])/theta[2])
  CvM(z)
}
