# Test functions: cdf.gamma, cdf.normal, cdf.weibull, cdf.logistic
# Input checking required for both x and parameters


#######################################################

# Test functions: AD.logistic.pvalue, AD.weibull.pvalue
# Test functions: CvM.logistic.pvalue, CvM.weibull.pvalue, CvM.poisson.pvalue.boostrap
# Function missing:

AD.logistic.pvalue(0.5) #play with eps and limit; statistics from 0.1 to 10
AD.weibull.pvalue(0.5)
AD.normal.pvalue(0.5) #e-7 very slow
AD.gamma.pvalue(0.5)
AD.laplace.pvalue(0.5) #e-7
AD.exp.pvalue(0.5)

CvM.logistic.pvalue(0.05)
CvM.weibull.pvalue(0.05)
CvM.normal.pvalue(0.05) #e-7
CvM.gamma.pvalue(0.05)
CvM.laplace.pvalue(0.05) #e-7
CvM.exp.pvalue(0.05)


CvM.poisson.pvalue.bootstrap(rpois(100,2)) #working

# Monte Carlo

check.AD.logistic = function(par,n,m){
  out <- matrix(0,m,1)
  pval <- matrix(0,m,1)
  for(i in 1:m){
    x=rlogis(n,par[1],par[2])
    out[i,]=AD.logistic(x)
    pval[i,]=AD.logistic.pvalue(out[i,])$P
  }
  hist(out,main="AD statistics for Logistic random sample")
  abline(v=mean(out), col="red")
  hist(pval,main="p value of AD")
  abline(v=mean(pval), col="red")
}
system.time(check.AD.logistic(c(2,1),100,100))

check.AD.weibull = function(par,n,m){
  out <- matrix(0,m,1)
  pval <- matrix(0,m,1)
  for(i in 1:m){
    x=rweibull(n,par[1],par[2])
    out[i,]=AD.weibull(x)
    pval[i,]=AD.weibull.pvalue(out[i,])$P
  }
  hist(out,main="AD statistics for Weibull random sample")
  abline(v=mean(out), col="red")
  hist(pval,main="p value of AD")
  abline(v=mean(pval), col="red")
}
system.time(check.AD.weibull(c(2,1),500,100))

check.AD.normal = function(par,n,m){
  out <- matrix(0,m,1)
  pval <- matrix(0,m,1)
  for(i in 1:m){
    x=rnorm(n,par[1],par[2])
    out[i,]=AD.normal(x)
    pval[i,]=AD.normal.pvalue(out[i,])$P
  }
  hist(out,main="AD statistics for Normal random sample")
  abline(v=mean(out), col="red")
  hist(pval,main="p value of AD")
  abline(v=mean(pval), col="red")
}
system.time(check.AD.normal(c(2,1),500,10))



