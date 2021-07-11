#' Watson.laplace
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
Watson.laplace = function(x){
  pars <- estimate.laplace(x)
  z <- cdf.laplace((x-par[1])/pars[2])
  Watson(z)
}
