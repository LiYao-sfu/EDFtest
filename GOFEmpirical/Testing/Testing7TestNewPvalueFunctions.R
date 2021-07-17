
library(GOFEmpirical)

#
#  Gamma shape 2
#
set.seed(1943)
n=100
M=1000
shape=2
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

logistic.dat = matrix(rlogis(n*M,0,100),nrow=M)
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
shape=2
laplace.dat = matrix(LaplacesDemon::rlaplace(n*M,100,100),nrow=M)
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
shape=2
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
#  Gamma shape 10
#
set.seed(1943)
n=100
M=1000
shape=10
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
#  Gamma shape 0.5
#
set.seed(1943)
n=100
M=1000
shape=0.5
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
