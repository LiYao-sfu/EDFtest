#' gof.weibull.bootstrap
#'
#' @param x
#' @param M
#'
#' @return
#' @export
#'
#' @examples
"gof.weibull.bootstrap" <-
  function(x, M=10000){
    a2 <- AD.weibull(x)
    w2 <- CvM.weibull(x)
    n <- length(x)
    pars <- estimate.weibull(x)
    alpha <- pars[1]
    beta <- pars[2]
    dat <- rweibull(n*M,shape=alpha,scale=beta)
    dat <- matrix(dat,nrow=M)
    a2vals <- apply(dat,1,AD.weibull)
    w2vals <- apply(dat,1,CvM.weibull)
    a.pv <- length(a2vals[a2vals>a2])/M
    w.pv <- length(w2vals[w2vals>w2])/M
    a2text <- paste("A squared is ", as.character(round(a2,4)))
    a2text <- paste(a2text,"P-value is ", as.character(round(a.pv,4)),"\n")
    w2text <- paste("W squared is ", as.character(round(w2,4)))
    w2text <- paste(w2text,"P-value is ", as.character(round(w.pv,4)),"\n")
    cat(a2text)
    cat(w2text)
    invisible(list(Asq = a2,Asq.pvalue =a.pv, Wsq=w2, Wsq.pvalue = w2))
  }
