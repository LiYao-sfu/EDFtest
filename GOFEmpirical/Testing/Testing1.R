# Test functions: estimate.normal; estimate.gamma; estimate.weibull; estimate.logistic

library(GOFEmpirical)

# Normal
x1=rnorm(1000,7,12)
estimate.normal(x1)
x1=rnorm(1000,5,999)
estimate.normal(x1)

# Gamma
x2=rgamma(1000,7,12)
estimate.gamma(x2)
x2=rgamma(1000,0.01,999)
estimate.gamma(x2)

# Logistic
x3=rlogis(1000,7,12)
estimate.logistic(x3)
x3=rlogis(1000,0.01,999)
estimate.logistic(x3)

# Weibull
x4=rweibull(100000,7,12)
estimate.weibull(x4)
x4=rweibull(100000,0.01,999)
estimate.weibull(x4)

# Laplace
#install.packages("ExtDist")
library(ExtDist)
x6=rLaplace(1000,0,1)
estimate.laplace(x6)
x6=rLaplace(1000,2,66)
estimate.laplace(x6)

# Exponential
x7=rexp(1000,1)
estimate.exp(x7)
x7=rexp(1000,100)
estimate.exp(x7)

par(mfrow=c(1,2))

check.estimation.normal = function(par,n,m){
  out <- matrix(0,m,2)
  for(i in 1:m){
  x=rnorm(n,par[1],par[2])
  out[i,]=estimate.normal(x)
  }
  hist(out[,1],main="Parameter 1")
  abline(v=par[1], col="red")
  hist(out[,2],main="Parameter 2")
  abline(v=par[2], col="red")
}
system.time(check.estimation.normal(c(7,13),100000,1000))

check.estimation.gamma = function(par,n,m){
  out <- matrix(0,m,2)
  for(i in 1:m){
    x=rgamma(n,par[1],par[2])
    out[i,]=estimate.gamma(x)
  }
  hist(out[,1],main=paste0("Shape: ",par[1]))
  abline(v=par[1], col="red")
  hist(1/out[,2],main=paste0("Scale: ",1/par[2]))
  abline(v=par[2], col="red")
}
system.time(check.estimation.gamma(c(7,13),100000,1000))#shape 0.5-50,


check.estimation.logistic = function(par,n,m){
  out <- matrix(0,m,2)
  for(i in 1:m){
    x=rlogis(n,par[1],par[2])
    out[i,]=estimate.logistic(x)
  }
  hist(out[,1],main="Parameter 1")
  abline(v=par[1], col="red")
  hist(out[,2],main="Parameter 2")
  abline(v=par[2], col="red")
}
system.time(check.estimation.logistic(c(7,13),100000,1000))

check.estimation.weibull = function(par,n,m){
  out <- matrix(0,m,2)
  for(i in 1:m){
    x=rweibull(n,par[1],par[2])
    out[i,]=estimate.weibull(x)
  }
  hist(out[,1],main=paste0("Shape: ",par[1]))
  abline(v=par[1], col="red")
  hist(out[,2],main=paste0("Scale: ",par[2]))
  abline(v=par[2], col="red")
}
system.time(check.estimation.weibull(c(7,13),100000,1000)) #set shape 0.5-50

check.estimation.laplace = function(par,n,m){
  out <- matrix(0,m,2)
  for(i in 1:m){
    x=rLaplace(n,par[1],par[2])
    out[i,]=estimate.laplace(x)
  }
  hist(out[,1],main=paste0("Location: ",par[1]))
  abline(v=par[1], col="red")
  hist(out[,2],main=paste0("Scale: ",par[2]))
  abline(v=par[2], col="red")
}
system.time(check.estimation.laplace(c(0,1),1000,1000))

check.estimation.exp = function(par,n,m){
  out <- matrix(0,m,1)
  for(i in 1:m){
    x=rexp(n,par[1])
    out[i,]=estimate.exp(x)
  }
  hist(out[,1],main=paste0("Scale: ",par[1]))
  abline(v=par[1], col="red")
}
system.time(check.estimation.exp(c(1),1000,1000))

##############################################################################

# Gamma and Weibull

par(mfrow=c(1,2))

for(i in 1:5){
  check.estimation.gamma(c(0.5*i,1),1000,1000)
}
check.estimation.gamma(c(5,1),1000,1000)
check.estimation.gamma(c(10,1),1000,1000)
check.estimation.gamma(c(50,1),1000,1000)


for(i in 1:5){
  check.estimation.weibull(c(0.5*i,1),1000,1000)
}
check.estimation.weibull(c(5,1),1000,1000)
check.estimation.weibull(c(10,1),1000,1000)
check.estimation.weibull(c(50,1),1000,1000)


##############################################################################







