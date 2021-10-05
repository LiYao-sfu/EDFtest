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
#' @param m
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
#' gof.gamma.bootstrap(x2,M=1000)
#'
#' x3=rlogis(n=1000,location=0,scale=1)
#' gof.logistic(x3)
#' gof.logistic.bootstrap(x3,M=1000)
#'
#' x4= rmutil::rlaplace(n=1000,m=0,s=1) 
#' gof.laplace(x4)
#' gof.laplace.bootstrap(x4,M=1000)
#' 
#' x5 = rweibull(100,1)
#' gof.weibull(x5)
#' gof.weibull(x5,print=TRUE,verbose=TRUE)
#' 
#' x6=rweibull(1000,1)
#' gof.weibull(x6)
#' gof.weibull.bootstrap(x6,M=1000)
#' 
#' x7 = rexp(100,1/2)
#' gof.exp(x7)
#' gof.exp(x7,print=TRUE,verbose=TRUE)
#' 
#' x8=rexp(1000,1/2)
#' gof.exp.bootstrap(x8,M=1000)
#' gof.exp(x8)

library(rmutil)

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
#' @rdname gof.normal.bootstrap
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
#' @rdname gof.gamma
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
#' @rdname gof.gamma.bootstrap
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
#' @rdname gof.logistic
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
  invisible(list(Asq = a2,Asq.pvalue =a.pv, Wsq=w2, Wsq.pvalue = w2,Usq=u2,Usq.pvalue=u2))
}
#' @export
#' @rdname gof.logistic.bootstrap

gof.laplace=function(x,print=TRUE,verbose=FALSE){
  
  pars=estimate.laplace(x)
  if(verbose){cat("Laplace parameter estimates", pars, "\n")}
  pit=rmutil::plaplace(x,m=pars[1],s=pars[2])
  if(verbose){cat("PITs are done \n \n")}
  w = CvM(pit)
  a = AD(pit)
  u = Watson(pit)
  w.p=CvM.laplace.pvalue(w)$P
  if(verbose){
    cat("Cramer von Mises P value output \n")
    print(w.p)
    cat("\n\n")
  }
  w.p=w.p$P
  a.p=AD.laplace.pvalue(a)$P
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
#' @rdname gof.laplace

gof.laplace.bootstrap<-function(x, M=10000){
  a2 <- AD.laplace(x)
  w2 <- CvM.laplace(x)
  u2 <- Watson.laplace(x)
  n <- length(x)
  pars <- estimate.laplace(x)
  med <- pars[1]
  MAD <- pars[2]
  dat <- rlaplace(n*M,location=med,scale=MAD)
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
  invisible(list(Asq = a2,Asq.pvalue =a.pv, Wsq=w2, Wsq.pvalue = w2, Usq=u2, Usq.pvalue = u2))
}

#' @export
#' @rdname gof.laplace.bootstrap

gof.weibull=function(x,print=TRUE,verbose=FALSE){
  
  pars=estimate.weibull(x)
  if(verbose){cat("Weibull parameter estimates", pars, "\n")}
  
  pit=pweibull(x,shape=pars[1],scale=pars[2])
  if(verbose){cat("PITs are done \n \n")}
  w = CvM(pit)
  a = AD(pit)
  u = Watson(pit)
  w.p=CvM.weibull.pvalue(w)$P
  if(verbose){
    cat("Cramer von Mises P value output \n")
    print(w.p)
    cat("\n\n")
  }
  w.p=w.p$P
  a.p=AD.weibull.pvalue(a)$P
  if(verbose){
    cat("Anderson-Darling P value output \n")
    print(a.p)
    cat("\n\n")
  }
  a.p=a.p$P
  u.p=Watson.weibull.pvalue(u)$P
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
#' @rdname gof.weibull

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
  invisible(list(Asq = a2,Asq.pvalue =a.pv, Wsq=w2, Wsq.pvalue = w2, Usq=u2, Usq.pvalue = u2))
}
#' @export
#' @rdname gof.weibull.bootstrap

gof.exp=function(x,print=TRUE,verbose=FALSE){
  pars=estimate.exp(x)
  if(verbose){cat("Exponential parameter estimates", pars, "\n")}

  pit=pexp(x,rate=1/pars)
  if(verbose){cat("PITs are done \n \n")}
  w = CvM(pit)
  a = AD(pit)
  u = Watson(pit)
  
  w.p=CvM.exp.pvalue(w)$P
  if(verbose){
    cat("Cramer von Mises P value output \n")
    print(w.p)
    cat("\n\n")
  }
  a.p=AD.exp.pvalue(a)$P
  if(verbose){
    cat("Anderson-Darling P value output \n")
    print(a.p)
    cat("\n\n")
  }
  u.p=Watson.exp.pvalue(a)$P
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
#' @export
#' @rdname gof.exp

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
  invisible(list(Asq = a2,Asq.pvalue =a.pv, Wsq=w2, Wsq.pvalue = w2, Usq=u2, Usq.pvalue = u2))
}
#' @export
#' @rdname gof.exp.bootstrap

gof.uniform.bootstrap=function(x,stat=y,m=10000){
  Z=sort(x)
  n=length(Z)
  r=matrix(runif(n*m),nrow=m)
  r = t(apply(r,1,sort))
  print(dim(r))
  if(stat==1)
  {
    Wsq=CvM(Z)
    Wsq_r=apply(r,1,CvM)
    pval=sum(Wsq_r>=Wsq)/m
  }
  else if(stat==2)
  {
    Usq=usq(Z)
    Usq_r=apply(r,1,usq)
    pval=sum(Usq_r>=Usq)/m
  }
  else
  {
    Asq=AD(Z)
    Asq_r=apply(r,1,AD)
    pval=sum(Asq_r>=Asq)/m
  }
  return(pval)
}
#' @export
#' @rdname gof.uniform.bootstrap
