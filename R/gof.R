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
#' @param stat
#'
#' @return Cramer-von Mises, Anderson-Darling and Watson's statistics and their P-values.
#' @export
#'
#' @examples
#' x0=runif(n=100,min=-1,max=1)
#' gof.uniform(x0)
#' gof.uniform.bootstrap(x0,M=1000)
#'
#' x1=rnorm(n=100,mean=0,sd=1)
#' gof.normal(x1)
#' gof.normal.bootstrap(x1,M=1000)
#'
#' x2=rgamma(n=100,shape=1,scale=1)
#' gof.gamma(x2)
#' gof.gamma.bootstrap(x2,M=1000)
#'
#' x3=rlogis(n=100,location=0,scale=1)
#' gof.logistic(x3)
#' gof.logistic.bootstrap(x3,M=1000)
#'
#' x4= rmutil::rlaplace(n=100,m=0,s=1)
#' gof.laplace(x4)
#' gof.laplace.bootstrap(x4,M=1000)
#'
#' x5=rweibull(n=100,shape=1,scale=1)
#' gof.weibull(x5)
#' gof.weibull.bootstrap(x5,M=1000)
#'
#' x6=rexp(n=100,rate=1/2)
#' gof.exp(x6)
#' gof.exp.bootstrap(x6,M=1000)
gof.uniform=function(x,print=TRUE,verbose=FALSE){
  #  Estimate the parameters
  pars=estimate.uniform(x)
  if(verbose){cat("Uniform parameter estimates", pars, "\n")}

  #  Compute the pit
  pit=punif(x,min=pars[1],max=pars[2])
  if(verbose){cat("PITs are done \n \n")}

  #  Compute two gof statistics
  w = CvM(pit)
  a = AD(pit)
  u = Watson(pit)

  #  Compute their p-values
  w.p=CvM.uniform.pvalue(w,verbose=verbose)$P
  if(verbose){
    cat("Cramer-von Mises P-value output \n")
    print(w.p)
    cat("\n\n")
  }
  a.p=AD.uniform.pvalue(a,verbose=verbose)$P
  if(verbose){
    cat("Anderson-Darling P-value output \n")
    print(a.p)
    cat("\n\n")
  }
  u.p=Watson.uniform.pvalue(u,verbose=verbose)$P
  if(verbose){
    cat("Watson P-value output \n")
    print(u.p)
    cat("\n\n")
  }
  if(print){
    cat("Cramer-von Mises statistic is ",w,"with P-value is ",w.p,"\n")
    cat("Anderson-Darling statistic is ",a,"with P-value is ",a.p,"\n")
    cat("Watson statistic is ",u,"with P-value is ",u.p,"\n")
  }
  invisible(list(Wsq=w,Wsq.pvalue=w.p,Asq=a,Asq.pvalue=a.p,Usq=u,Usq.pvalue=u.p))
}

#' @export
#' @rdname gof.uniform
gof.uniform.bootstrap<-function(x, M=10000){
  a2 <- AD.uniform(x)
  w2 <- CvM.uniform(x)
  u2 <- Watson.uniform(x)
  n <- length(x)
  pars <- estimate.uniform(x)
  mi <- pars[1]
  ma <- pars[2]
  dat <- runif(n*M,min=mi,max=ma)
  dat <- matrix(dat,nrow=M)
  a2vals <- apply(dat,1,AD.uniform)
  w2vals <- apply(dat,1,CvM.uniform)
  u2vals <- apply(dat,1,Watson.uniform)
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
  invisible(list(Wsq=w2,Wsq.pvalue=w.pv,Asq=a2,Asq.pvalue=a.pv,Usq=u2,Usq.pvalue=u.pv))
}

#' @export
#' @rdname gof.uniform
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
  w.p=CvM.normal.pvalue(w,verbose=verbose)$P
  if(verbose){
    cat("Cramer-von Mises P-value output \n")
    print(w.p)
    cat("\n\n")
  }
  a.p=AD.normal.pvalue(a,verbose=verbose)$P
  if(verbose){
    cat("Anderson-Darling P-value output \n")
    print(a.p)
    cat("\n\n")
  }
  u.p=Watson.normal.pvalue(u,verbose=verbose)$P
  if(verbose){
    cat("Watson P-value output \n")
    print(u.p)
    cat("\n\n")
  }
  if(print){
    cat("Cramer-von Mises statistic is ",w,"with P-value is ",w.p,"\n")
    cat("Anderson-Darling statistic is ",a,"with P-value is ",a.p,"\n")
    cat("Watson statistic is ",u,"with P-value is ",u.p,"\n")
  }
  invisible(list(Wsq=w,Wsq.pvalue=w.p,Asq=a,Asq.pvalue=a.p,Usq=u,Usq.pvalue=u.p))
}

#' @export
#' @rdname gof.uniform
gof.normal.bootstrap<-function(x, M=10000){
  pars <- estimate.normal(x)
  u <- cdf.normal(x,parameter=pars)
  a2 <- AD(u)
  w2 <- CvM(u)
  u2 <- Watson(u)
  n <- length(x)
  xbar <- pars[1]
  s <- pars[2]
  dat <- rnorm(n*M,mean=xbar,sd=s)
  dat <- matrix(dat,nrow=M)
  uvales <- apply(dat,1,function(x){cdf.normal(x,parameter=estimate.normal(x))})
  a2vals <- apply(uvales,1,AD)
  w2vals <- apply(uvales,1,CvM)
  u2vals <- apply(uvales,1,Watson)
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
  invisible(list(Wsq=w2,Wsq.pvalue=w.pv,Asq=a2,Asq.pvalue=a.pv,Usq=u2,Usq.pvalue=u.pv))
}

#' @export
#' @rdname gof.uniform
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
    cat("Cramer-von Mises P-value output \n")
    print(w.p)
    cat("\n\n")
  }
  a.p=AD.gamma.pvalue(a,shape=pars[1])$P
  if(verbose){
    cat("Anderson-Darling P-value output \n")
    print(a.p)
    cat("\n\n")
  }
  u.p=Watson.gamma.pvalue(u,shape=pars[1])$P
  if(verbose){
    cat("Anderson-Darling P-value output \n")
    print(u.p)
    cat("\n\n")
  }
  if(print){
    cat("Cramer-von Mises statistic is ",w,"with P-value is ",w.p,"\n")
    cat("Anderson-Darling statistic is ",a,"with P-value is ",a.p,"\n")
    cat("Watson statistic is ",u,"with P-value is ",u.p,"\n")
  }
  invisible(list(Wsq=w,Wsq.pvalue=w.p,Asq=a,Asq.pvalue=a.p,Usq=u,Usq.pvalue=u.p))
}

#' @export
#' @rdname gof.uniform
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
  invisible(list(Wsq=w2,Wsq.pvalue=w.pv,Asq=a2,Asq.pvalue=a.pv,Usq=u2,Usq.pvalue=u.pv))
}

#' @export
#' @rdname gof.uniform
gof.logistic=function(x,print=TRUE,verbose=FALSE){
  pars=estimate.logistic(x)
  if(verbose){cat("log-Logistic parameter estimates", pars, "\n")}
  pit=plogis(x,location=pars[1],scale=pars[2])
  if(verbose){cat("PITs are done \n \n")}
  w = CvM(pit)
  a = AD(pit)
  u = Watson(pit)
  if(verbose){cat("Statistics are ",w,a,u," \n \n")}
  w.p=CvM.logistic.pvalue(w,verbose=verbose)$P
  if(verbose){
    cat("Cramer-von Mises P-value output \n")
    print(w.p)
    cat("\n\n")
  }
  a.p=AD.logistic.pvalue(a,verbose=verbose)$P
  if(verbose){
    cat("Anderson-Darling P-value output \n")
    print(a.p)
    cat("\n\n")
  }
  u.p=Watson.logistic.pvalue(u,verbose=verbose)$P
  if(verbose){
    cat("Watson P-value output \n")
    print(u.p)
    cat("\n\n")
  }
  if(print){
    cat("Cramer-von Mises statistic is ",w,"with P-value is ",w.p,"\n")
    cat("Anderson-Darling statistic is ",a,"with P-value is ",a.p,"\n")
    cat("Watson statistic is ",u,"with P-value is ",u.p,"\n")
  }
  invisible(list(Wsq=w,Wsq.pvalue=w.p,Asq=a,Asq.pvalue=a.p,Usq=u,Usq.pvalue=u.p))
}

#' @export
#' @rdname gof.uniform
gof.logistic.bootstrap=function(x, M=10000){
  a2 <- AD.logistic(x)
  w2 <- CvM.logistic(x)
  u2 <- Watson.logistic(x)
  n <- length(x)
  pars <- estimate.logistic(x)
  alpha <- pars[1]
  beta <- pars[2]
  dat <- rlogis(n*M,location=alpha,scale=beta)
  dat <- matrix(dat,nrow=M)
  a2vals <- apply(dat,1,AD.logistic)
  w2vals <- apply(dat,1,CvM.logistic)
  u2vals <- apply(dat,1,Watson.logistic)
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
  invisible(list(Wsq=w2,Wsq.pvalue=w.pv,Asq=a2,Asq.pvalue=a.pv,Usq=u2,Usq.pvalue=u.pv))
}

#' @export
#' @rdname gof.uniform
gof.laplace=function(x,print=TRUE,verbose=FALSE){
  pars=estimate.laplace(x,use.sd=FALSE)
  if(verbose){cat("Laplace parameter estimates", pars, "\n")}
  pit=cdf.laplace(x,theta =pars)
  if(verbose){cat("PITs are done \n \n")}
  w = CvM(pit)
  a = AD(pit)
  u = Watson(pit)
  w.p=CvM.laplace.pvalue(w,verbose=verbose)$P
  if(verbose){
    cat("Cramer-von Mises P-value output \n")
    print(w.p)
    cat("\n\n")
  }
  a.p=AD.laplace.pvalue(a,verbose=verbose)$P
  if(verbose){
    cat("Anderson-Darling P-value output \n")
    print(a.p)
    cat("\n\n")
  }
  u.p=Watson.logistic.pvalue(u,verbose=verbose)$P
  if(verbose){
    cat("Watson P-value output \n")
    print(u.p)
    cat("\n\n")
  }
  if(print){
    cat("Cramer-von Mises statistic is ",w,"with P-value is ",w.p,"\n")
    cat("Anderson-Darling statistic is ",a,"with P-value is ",a.p,"\n")
    cat("Watson statistic is ",u,"with P-value is ",u.p,"\n")
  }
  invisible(list(Wsq=w,Wsq.pvalue=w.p,Asq=a,Asq.pvalue=a.p,Usq=u,Usq.pvalue=u.p))
}

#' @export
#' @rdname gof.uniform
gof.laplace.bootstrap<-function(x, M=10000){
  a2 <- AD.laplace(x)
  w2 <- CvM.laplace(x)
  u2 <- Watson.laplace(x)
  n <- length(x)
  pars <- estimate.laplace(x)
  med <- pars[1]
  MAD <- pars[2]
  dat <- rmutil::rlaplace(n*M,m=med,s=MAD)
  dat <- matrix(dat,nrow=M)
  a2vals <- apply(dat,1,AD.laplace)
  w2vals <- apply(dat,1,CvM.laplace)
  u2vals <- apply(dat,1,Watson.laplace)
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
  invisible(list(Wsq=w2,Wsq.pvalue=w.pv,Asq=a2,Asq.pvalue=a.pv,Usq=u2,Usq.pvalue=u.pv))
}

#' @export
#' @rdname gof.uniform
gof.weibull=function(x,print=TRUE,verbose=FALSE){
  pars=estimate.weibull(x)
  if(verbose){cat("Weibull parameter estimates", pars, "\n")}
  pit=pweibull(x,shape=pars[1],scale=pars[2])
  if(verbose){cat("PITs are done \n \n")}
  w = CvM(pit)
  a = AD(pit)
  u = Watson(pit)
  w.p=CvM.weibull.pvalue(w,verbose=verbose)$P
  if(verbose){
    cat("Cramer-von Mises P-value output \n")
    print(w.p)
    cat("\n\n")
  }
  a.p=AD.weibull.pvalue(a,verbose=verbose)$P
  if(verbose){
    cat("Anderson-Darling P-value output \n")
    print(a.p)
    cat("\n\n")
  }
  u.p=Watson.weibull.pvalue(u,verbose=verbose)$P
  if(verbose){
    cat("Watson P-value output \n")
    print(u.p)
    cat("\n\n")
  }
  if(print){
    cat("Cramer-von Mises statistic is ",w,"with P-value is ",w.p,"\n")
    cat("Anderson-Darling statistic is ",a,"with P-value is ",a.p,"\n")
    cat("Watson statistic is ",u,"with P-value is ",u.p,"\n")
  }
  invisible(list(Wsq=w,Wsq.pvalue=w.p,Asq=a,Asq.pvalue=a.p,Usq=u,Usq.pvalue=u.p))
}

#' @export
#' @rdname gof.uniform
gof.weibull.bootstrap<-function(x, M=10000){
  a2 <- AD.weibull(x)
  w2 <- CvM.weibull(x)
  u2 <- Watson.weibull(x)
  n <- length(x)
  pars <- estimate.weibull(x)
  alpha <- pars[1]
  beta <- pars[2]
  dat <- rweibull(n*M,shape=alpha,scale=beta)
  dat <- matrix(dat,nrow=M)
  a2vals <- apply(dat,1,AD.weibull)
  w2vals <- apply(dat,1,CvM.weibull)
  u2vals <- apply(dat,1,Watson.weibull)
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
  invisible(list(Wsq=w2,Wsq.pvalue=w.pv,Asq=a2,Asq.pvalue=a.pv,Usq=u2,Usq.pvalue=u.pv))
}

#' @export
#' @rdname gof.uniform
gof.exp=function(x,print=TRUE,verbose=FALSE){
  pars=estimate.exp(x)
  if(verbose){cat("Exponential parameter estimates", pars, "\n")}
  pit=pexp(x,rate=1/pars)
  if(verbose){cat("PITs are done \n \n")}
  w = CvM(pit)
  a = AD(pit)
  u = Watson(pit)
  w.p=CvM.exp.pvalue(w,verbose=verbose)$P
  if(verbose){
    cat("Cramer-von Mises P-value output \n")
    print(w.p)
    cat("\n\n")
  }
  a.p=AD.exp.pvalue(a,verbose=verbose)$P
  if(verbose){
    cat("Anderson-Darling P-value output \n")
    print(a.p)
    cat("\n\n")
  }
  u.p=Watson.exp.pvalue(u,verbose=verbose)$P
  if(verbose){
    cat("Watson P-value output \n")
    print(u.p)
    cat("\n\n")
  }
  if(print){
    cat("Cramer-von Mises statistic is ",w,"with P-value is ",w.p,"\n")
    cat("Anderson-Darling statistic is ",a,"with P-value is ",a.p,"\n")
    cat("Watson statistic is ",u,"with P-value is ",u.p,"\n")
  }
  invisible(list(Wsq=w,Wsq.pvalue=w.p,Asq=a,Asq.pvalue=a.p,Usq=u,Usq.pvalue=u.p))
}

#' @export
#' @rdname gof.uniform
gof.exp.bootstrap<-function(x, M=10000){
  a2 <- AD.exp(x)
  w2 <- CvM.exp(x)
  u2 <- Watson.exp(x)
  n <- length(x)
  pars <- estimate.exp(x)
  scale <- pars
  dat <- rexp(n*M,rate=1/scale)
  dat <- matrix(dat,nrow=M)
  a2vals <- apply(dat,1,AD.exp)
  w2vals <- apply(dat,1,CvM.exp)
  u2vals <- apply(dat,1,Watson.exp)
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
  invisible(list(Wsq=w2,Wsq.pvalue=w.pv,Asq=a2,Asq.pvalue=a.pv,Usq=u2,Usq.pvalue=u.pv))
}

