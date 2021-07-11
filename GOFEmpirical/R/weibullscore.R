#' weibullscore
#'
#' @param a
#' @param x
#'
#' @return
#' @export
#'
#' @examples
"weibullscore" <-
  function(a,x){
    ml <- mean(log(x))
    pow <- t(outer(x,a,"^"))
    one <- rep(1,length)
    mla <- pow%*%log(x)
    ma <- pow%*%one
    ml - mla/ma+log(a) -log(ma)/a +1/a
  }
