#' EDF Goodness-of-Fit Tests for a Given Distribution
#'
#' @description
#' This function takes in an i.i.d. random sample, use MLE to estimate Normal
#' parameters, compute probability integral transforms, and computes Cramer-von Mises,
#' Anderson-Darling and Watson's statistics and their P-values using either `imhof`
#' method or bootstrap.
#'
#' @details
#'
#' @param x A random sample.
#' @param print Logical; if TRUE print both statistics and P-values.
#' @param verbose verbose Logical; if TRUE, print warning messages.
#' @param M Number of bootstrap, 10000 by default.
#'
#' @return Cramer-von Mises, Anderson-Darling and Watson's statistics and their P-values.
#' @export
#'
#' @examples
#' x1=rnorm(n=1000,mean=0,sd=1)
#' gof.normal(x1)
#' gof.normal.bootstrap(x1,M=1000)
#'
#' x2=rgamma(n=1000,shape=1,scale=1)
#' gof.gamma(x)
#' gof.normal.bootstrap(x,M=1000)
#'
#' x3=rlogis(n=1000,location=0,scale=1)
#' usq3 = Watson.logistic(x)
#' Watson.logistic.pvalue(usq3)
#'
#' x4= rmutil::rlaplace(n=1000,m=0,s=1)
#' gof.laplace(x4)
#' gof.laplace.bootstrap(x4,M=1000)
gof.normal=function(x,print=TRUE,verbose=FALSE){
  #  Estimate the parameters
  pars=estimate.normal(x)
  if(verbose){cat("Normal parameter estimates", pars, "\n")}

  #  Compute the pit
  pit=pnorm(x,mean=pars[1],sd=pars[2])
  if(verbose){cat("PITs are done \n \n")}

  #  Compute two gof statistics
  w = CvM(pit)
  a = AD(pit)
  u = Watson(pit)

  #  Compute their p-values
  w.p=CvM.normal.pvalue(w)$P
  if(verbose){
    cat("Cramer von Mises P value output \n")
    print(w.p)
    cat("\n\n")
  }
  a.p=AD.normal.pvalue(a)$P
  if(verbose){
    cat("Anderson-Darling P value output \n")
    print(a.p)
    cat("\n\n")
  }
  u.p=Watson.normal.pvalue(a)$P
  if(verbose){
    cat("Watson P value output \n")
    print(u.p)
    cat("\n\n")
  }
  if(print){
    cat("Cramer-von Mises statistic is ",w,"with P-value is ",w.p,"\n")
    cat("Anderson-Darling statistic is ",a,"with P-value is ",a.p,"\n")
    cat("Watson statistic is ",u,"with P-value is ",u.p,"\n")
  }
  invisible(list(w=w,w.p=w.p,a=a,a.p=a.p,u=u,u.p=u.p))
}

# Make bootstrap functions efficient!

#' @export
#' @rdname gof.normal
gof.normal.bootstrap<-function(x, M=10000){
  a2 <- AD.normal(x)
  w2 <- CvM.normal(x)
  u2 <- Watson.normal(x)
  n <- length(x)
  pars <- estimate.normal(x)
  xbar <- pars[1]
  s <- pars[2]
  dat <- rnorm(n*M,mean=xbar,sd=s)
  dat <- matrix(dat,nrow=M)
  a2vals <- apply(dat,1,AD.normal)
  w2vals <- apply(dat,1,CvM.normal)
  u2vals <- apply(dat,1,Watson.normal)
  a.pv <- length(a2vals[a2vals>a2])/M
  w.pv <- length(w2vals[w2vals>w2])/M
  u.pv <- length(u2vals[u2vals>u2])/M
  w2text <- paste("Cramer-von Mises statistic is ", as.character(round(w2,7)))
  w2text <- paste(w2text,"with P-value is ", as.character(round(w.pv,7)),"\n")
  a2text <- paste("Anderson-Darling statistic is ", as.character(round(a2,7)))
  a2text <- paste(a2text,"with P-value is ", as.character(round(a.pv,7)),"\n")
  u2text <- paste("Watson statistic is ", as.character(round(u2,7)))
  u2text <- paste(u2text,"with P-value is ", as.character(round(u.pv,7)),"\n")
  cat(w2text)
  cat(a2text)
  cat(u2text)
  invisible(list(Asq = a2,Asq.pvalue =a.pv, Wsq=w2, Wsq.pvalue = w2, Usq=u2, Usq.pvalue = u2))
}

#' @export
#' @rdname gof.normal
gof.gamma=function(x,print=TRUE,verbose=FALSE){
  pars=estimate.gamma(x)
  if(verbose){cat("Gamma parameter estimates", pars, "\n")}
  pit=pgamma(x,shape=pars[1],scale=pars[2])
  if(verbose){cat("PITs are done \n \n")}
  w = CvM(pit)
  a = AD(pit)
  u = Watson(pit)
  w.p=CvM.gamma.pvalue(w,shape=pars[1])$P
  if(verbose){
    cat("Cramer von Mises P value output \n")
    print(w.p)
    cat("\n\n")
  }
  a.p=AD.gamma.pvalue(a,shape=pars[1])$P
  if(verbose){
    cat("Anderson-Darling P value output \n")
    print(a.p)
    cat("\n\n")
  }
  u.p=AD.gamma.pvalue(u,shape=pars[1])$P
  if(verbose){
    cat("Anderson-Darling P value output \n")
    print(u.p)
    cat("\n\n")
  }
  if(print){
    cat("Cramer-von Mises statistic is ",w,"with P-value is ",w.p,"\n")
    cat("Anderson-Darling statistic is ",a,"with P-value is ",a.p,"\n")
    cat("Watson statistic is ",u,"with P-value is ",u.p,"\n")
  }
  invisible(list(w=w,w.p=w.p,a=a,a.p=a.p,u=u,u.p=u.p))
}

#' @export
#' @rdname gof.normal
gof.gamma.bootstrap<-function(x, M=10000){
  a2 <- AD.gamma(x)
  w2 <- CvM.gamma(x)
  u2 <- Watson.gamma(x)
  n <- length(x)
  pars <- estimate.gamma(x)
  alpha <- pars[1]
  beta <- pars[2]
  dat <- rgamma(n*M,shape=alpha,scale=beta)
  dat <- matrix(dat,nrow=M)
  a2vals <- apply(dat,1,AD.gamma)
  w2vals <- apply(dat,1,CvM.gamma)
  u2vals <- apply(dat,1,Watson.gamma)
  a.pv <- length(a2vals[a2vals>a2])/M
  w.pv <- length(w2vals[w2vals>w2])/M
  u.pv <- length(u2vals[w2vals>u2])/M
  w2text <- paste("Cramer-von Mises statistic is ", as.character(round(w2,7)))
  w2text <- paste(w2text,"with P-value is ", as.character(round(w.pv,7)),"\n")
  a2text <- paste("Anderson-Darling statistic is ", as.character(round(a2,7)))
  a2text <- paste(a2text,"with P-value is ", as.character(round(a.pv,7)),"\n")
  u2text <- paste("Watson statistic is ", as.character(round(u2,7)))
  u2text <- paste(u2text,"with P-value is ", as.character(round(u.pv,7)),"\n")
  cat(w2text)
  cat(a2text)
  cat(u2text)
  invisible(list(Asq = a2,Asq.pvalue =a.pv, Wsq=w2, Wsq.pvalue = w2, Usq=u2, Usq.pvalue = u2))
}

#' @export
#' @rdname gof.normal
gof.logistic=function(x,print=TRUE,verbose=FALSE){
  pars=estimate.logistic(x)
  if(verbose){cat("log-Logistic parameter estimates", pars, "\n")}
  pit=plogis(x,location=pars[1],scale=pars[2])
  if(verbose){cat("PITs are done \n \n")}
  w = CvM(pit)
  a = AD(pit)
  u = Watson(pit)
  if(verbose){cat("Statistics are ",w,a,u," \n \n")}
  w.p=CvM.logistic.pvalue(w,verbose=verbose)
  if(verbose){
    cat("Cramer von Mises P value output \n")
    print(w.p)
    cat("\n\n")
  }
  w.p=w.p$P
  a.p=AD.logistic.pvalue(a)$P
  if(verbose){
    cat("Anderson-Darling P value output \n")
    print(a.p)
    cat("\n\n")
  }
  a.p=a.p$P
  u.p=Watson.logistic.pvalue(u)$P
  if(verbose){
    cat("Watson P value output \n")
    print(u.p)
    cat("\n\n")
  }
  u.p=u.p$P
  if(print){
    cat("Cramer-von Mises statistic is ",w,"with P-value is ",w.p,"\n")
    cat("Anderson-Darling statistic is ",a,"with P-value is ",a.p,"\n")
    cat("Watson statistic is ",u,"with P-value is ",u.p,"\n")
  }
  invisible(list(w=w,w.p=w.p,a=a,a.p=a.p,u=u,u.p=u.p))
}

#' @export
#' @rdname gof.normal
gof.logistic.bootstrap=function(x, M=10000){
  a2 <- AD.logistic(x)
  w2 <- CvM.logistic(x)
  n <- length(x)
  pars <- estimate.logistic(x)
  alpha <- pars[1]
  beta <- pars[2]
  dat <- rlogis(n*M,location=alpha,scale=beta)
  dat <- matrix(dat,nrow=M)
  a2vals <- apply(dat,1,AD.logistic)
  w2vals <- apply(dat,1,CvM.logistic)
  a.pv <- length(a2vals[a2vals>a2])/M
  w.pv <- length(w2vals[w2vals>w2])/M
  w2text <- paste("Cramer-von Mises statistic is ", as.character(round(w2,7)))
  w2text <- paste(w2text,"with P-value is ", as.character(round(w.pv,7)),"\n")
  a2text <- paste("Anderson-Darling statistic is ", as.character(round(a2,7)))
  a2text <- paste(a2text,"with P-value is ", as.character(round(a.pv,7)),"\n")
  cat(w2text)
  cat(a2text)
  invisible(list(Asq = a2,Asq.pvalue =a.pv, Wsq=w2, Wsq.pvalue = w2))
}
