#' score.gamma3profile
#'
#' @param x
#' @param gamma
#'
#' @return
#' @export
#'
#' @examples
score.gamma3profile <- function(x,gamma){
  n <- length(x)
  M <- outer(x,gamma,"-")
  T1 <- apply(M,2,sum)
  T3 <- apply(1/M,2,sum)
  T2 <- apply(log(M),2, sum)
  alpha <- T1*T3/(n^2-T1*T3)
  beta <- T1/(n*alpha)
  U1 <- -n*log(beta)+T2-n*digamma(alpha)
  U1
}
