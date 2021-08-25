
library(GOFEmpirical)

#
#  Gamma
#
set.seed(1943)
n=100
M=1000

shape=1
gamma.dat = matrix(rgamma(n*M,shape=shape),nrow=M)
pars = apply(gamma.dat,1,estimate.gamma)
shape = pars[1,]
W.gamma = apply(gamma.dat,1,CvM.gamma)
A.gamma = apply(gamma.dat,1,AD.gamma)
hist(W.gamma,breaks=50)
hist(A.gamma,breaks=50)

PW.gamma = W.gamma
PA.gamma = A.gamma
for( i in 1:M){ #c(1:909,911:950,952:1000)
  PW.gamma[i] = CvM.gamma.pvalue(W.gamma[i],shape=shape[i])$P
  PA.gamma[i] = AD.gamma.pvalue(A.gamma[i],shape=shape[i])$P
}

hist(PW.gamma,breaks=20)
hist(PA.gamma,breaks=20)
plot(PW.gamma,PA.gamma, pch='.')
AD(PW.gamma)
AD(PA.gamma)

#
#  Logistic
#

set.seed(1943)
n=100
M=1000

shape=0.01

logistic.dat = matrix(rlogis(n*M,0,1),nrow=M)
pars = t(apply(logistic.dat,1,estimate.logistic))
W.logistic = apply(logistic.dat,1,CvM.logistic)
A.logistic = apply(logistic.dat,1,AD.logistic)
hist(W.logistic,breaks=50)
hist(A.logistic,breaks=50)

PW.logistic = W.logistic
PA.logistic = A.logistic
for( i in 6:M){ #c(1:909,911:950,952:1000)
  PW.logistic[i] = CvM.logistic.pvalue(W.logistic[i])$P
  PA.logistic[i] = AD.logistic.pvalue(A.logistic[i])$P
}

hist(PW.logistic,breaks=20)
hist(PA.logistic,breaks=20)
plot(PW.logistic,PA.logistic, pch='.')
AD(PW.logistic)
AD(PA.logistic)

#
#  Laplace
#

set.seed(1943)
n=100
M=1000
laplace.dat = matrix(LaplacesDemon::rlaplace(n*M,5,2),nrow=M)
pars = t(apply(laplace.dat,1,estimate.laplace))
W.laplace = apply(laplace.dat,1,CvM.laplace)
A.laplace = apply(laplace.dat,1,AD.laplace)
hist(W.laplace,breaks=50)
hist(A.laplace,breaks=50)

PW.laplace = W.laplace
PA.laplace = A.laplace
for( i in 1:M){ #c(1:909,911:950,952:1000)
  PW.laplace[i] = CvM.laplace.pvalue(W.laplace[i])$P
  PA.laplace[i] = AD.laplace.pvalue(A.laplace[i])$P
}

hist(PW.laplace,breaks=20)
hist(PA.laplace,breaks=20)
plot(PW.laplace,PA.laplace, pch='.')
AD(PW.laplace)
AD(PA.laplace)

#
#  Normal
#

set.seed(1943)
n=100
M=1000
normal.dat = matrix(rnorm(n*M),nrow=M)
pars = t(apply(normal.dat,1,estimate.normal))
W.normal = apply(normal.dat,1,CvM.normal)
A.normal = apply(normal.dat,1,AD.normal)
hist(W.normal,breaks=50)
hist(A.normal,breaks=50)

PW.normal = W.normal
PA.normal = A.normal
for( i in 1:M){ #c(1:909,911:950,952:1000)
  PW.normal[i] = CvM.normal.pvalue(W.normal[i])$P
  PA.normal[i] = AD.normal.pvalue(A.normal[i])$P
}

hist(PW.normal,breaks=20)
hist(PA.normal,breaks=20)
plot(PW.normal,PA.normal, pch='.')
AD(PW.normal)
AD(PA.normal)



#
#  Weibull
#
set.seed(1943)
n=100
M=1000
weibull.dat = matrix(rweibull(n*M,shape=50),nrow=M)
pars = apply(weibull.dat,1,estimate.weibull)
shape = pars[1,]
W.weibull = apply(weibull.dat,1,CvM.weibull)
A.weibull = apply(weibull.dat,1,AD.weibull)
#hist(W.weibull,breaks=50)
#hist(A.weibull,breaks=50)

PW.weibull = W.weibull
PA.weibull = A.weibull
for( i in 1:M){ #c(1:909,911:950,952:1000)
  PW.weibull[i] = CvM.weibull.pvalue(W.weibull[i])$P
  PA.weibull[i] = AD.weibull.pvalue(A.weibull[i])$P
}

hist(PW.weibull,breaks=20)
hist(PA.weibull,breaks=20)
plot(PW.weibull,PA.weibull, pch='.')
AD(PW.weibull)
AD(PA.weibull)



#
#  Exponential
#
set.seed(1943)
n=100
M=1000
exp.dat = matrix(rexp(n*M,rate=0.001),nrow=M)
pars = apply(exp.dat,1,estimate.exp)
rate = pars
W.exp = apply(exp.dat,1,CvM.exp)
A.exp = apply(exp.dat,1,AD.exp)
hist(W.exp,breaks=50)
hist(A.exp,breaks=50)

PW.exp = W.exp
PA.exp = A.exp
for( i in 1:M){ #c(1:909,911:950,952:1000)
  PW.exp[i] = CvM.exp.pvalue(W.exp[i])$P
  PA.exp[i] = AD.exp.pvalue(A.exp[i])$P
}

hist(PW.exp,breaks=20)
hist(PA.exp,breaks=20)
plot(PW.exp,PA.exp, pch='.')
AD(PW.exp)
AD(PA.exp)
