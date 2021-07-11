# Test functions: AD.gamma, AD.logistic, AD.normal, AD.weibull
# Adding functions: AD.exp, AD.laplace

# Missing function: AD.Poisson
# Answer: No AD for Poisson!!

check.AD.gamma = function(par,n,m){
  out <- matrix(0,m,1)
  for(i in 1:m){
    x=rgamma(n,par[1],par[2])
    out[i,]=AD.gamma(x)
  }
  hist(out,main="AD statistics for Gamma random sample")
  abline(v=mean(out), col="red")
}
system.time(check.AD.gamma(c(7,13),1000,100))

check.AD.logistic = function(par,n,m){
  out <- matrix(0,m,1)
  for(i in 1:m){
    x=rlogis(n,par[1],par[2])
    out[i,]=AD.logistic(x)
  }
  hist(out,main="AD statistics for Logistic random sample")
  abline(v=mean(out), col="red")
}
system.time(check.AD.logistic(c(7,13),1000,100))

check.AD.normal = function(par,n,m){
  out <- matrix(0,m,1)
  for(i in 1:m){
    x=rnorm(n,par[1],par[2])
    out[i,]=AD.normal(x)
  }
  hist(out,main="AD statistics for Normal random sample")
  abline(v=mean(out), col="red")
}
system.time(check.AD.normal(c(7,13),1000,100))

check.AD.weibull = function(par,n,m){
  out <- matrix(0,m,1)
  for(i in 1:m){
    x=rweibull(n,par[1],par[2])
    out[i,]=AD.weibull(x)
  }
  hist(out,main="AD statistics for Weibull random sample")
  abline(v=mean(out), col="red")
}
system.time(check.AD.weibull(c(7,13),1000,100))

check.AD.laplace = function(par,n,m){
  out <- matrix(0,m,1)
  for(i in 1:m){
    x=rLaplace(n,par[1],par[2])
    out[i,]=AD.laplace(x)
  }
  hist(out,main="AD statistics for Laplace random sample")
  abline(v=mean(out), col="red")
}
system.time(check.AD.laplace(c(0,1),1000,1000))

check.AD.exp = function(par,n,m){
  out <- matrix(0,m,1)
  for(i in 1:m){
    x=rexp(n,par[1])
    out[i,]=AD.exp(x)
  }
  hist(out,main="AD statistics for Exponential random sample")
  abline(v=mean(out), col="red")
}
system.time(check.AD.exp(c(1),1000,1000)) #Bug in AD.exp

###################################################################

# Test functions: CvM.gamma, CvM.logistic, CvM.normal, CvM.weibull, CvM.Poisson
# Adding functions: CvM.laplace, CvM.exponential

check.CvM.gamma = function(par,n,m){
  out <- matrix(0,m,1)
  for(i in 1:m){
    x=rgamma(n,par[1],par[2])
    out[i,]=CvM.gamma(x)
  }
  hist(out,main="CvM statistics for Gamma random sample")
  abline(v=mean(out), col="red")
}
system.time(check.CvM.gamma(c(7,13),1000,100))

check.CvM.logistic = function(par,n,m){
  out <- matrix(0,m,1)
  for(i in 1:m){
    x=rlogis(n,par[1],par[2])
    out[i,]=CvM.logistic(x)
  }
  hist(out,main="CvM statistics for Logistic random sample")
  abline(v=mean(out), col="red")
}
system.time(check.CvM.logistic(c(7,13),1000,100))

check.CvM.normal = function(par,n,m){
  out <- matrix(0,m,1)
  for(i in 1:m){
    x=rnorm(n,par[1],par[2])
    out[i,]=CvM.normal(x)
  }
  hist(out,main="CvM statistics for Normal random sample")
  abline(v=mean(out), col="red")
}
system.time(check.CvM.normal(c(7,13),1000,100))

check.CvM.weibull = function(par,n,m){
  out <- matrix(0,m,1)
  for(i in 1:m){
    x=rweibull(n,par[1],par[2])
    out[i,]=CvM.weibull(x)
  }
  hist(out,main="CvM statistics for Weibull random sample")
  abline(v=mean(out), col="red")
}
system.time(check.CvM.weibull(c(7,13),1000,100))

check.CvM.laplace = function(par,n,m){
  out <- matrix(0,m,1)
  for(i in 1:m){
    x=rLaplace(n,par[1],par[2])
    out[i,]=CvM.laplace(x)
  }
  hist(out,main="CvM statistics for Laplace random sample")
  abline(v=mean(out), col="red")
}
system.time(check.CvM.laplace(c(0,1),1000,1000))


check.CvM.exp = function(par,n,m){
  out <- matrix(0,m,1)
  for(i in 1:m){
    x=rexp(n,par[1])
    out[i,]=CvM.exp(x)
  }
  hist(out,main="CvM statistics for Exponential random sample")
  abline(v=mean(out), col="red")
}
system.time(check.CvM.exp(c(1),1000,1000))


# Add CvM.Poisson. Ask Richard: Do we need estimate.poisson? Answer: YES!

check.CvM.poisson = function(par,n,m){
  out <- matrix(0,m,1)
  for(i in 1:m){
    x=rpois(n,par)
    out[i,]=CvM.poisson(x)
  }
  hist(out,main="CvM statistics for Poisson random sample")
  abline(v=mean(out), col="red")
}
system.time(check.CvM.poisson(2,10000,1000))
system.time(check.CvM.poisson(0.5,10000,1000))
system.time(check.CvM.poisson(10,10000,1000))





