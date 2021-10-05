cdf.normal.user = function(x,theta){
  pnorm(x,mean=theta[1],sd=theta[2])
}

score.normal.user = function(x,theta){
  sig=theta[2]
  mu=theta[1]
  s.mean= (x-mu)/sig
  s.sd= s.mean^2/sig-length(x)/sig
  cbind(s.mean/sig,s.sd)
}
M=200
n=50
pval <- matrix(0,M,6)
colnames(pval) <- c("CvM","AD","Watson","CvM_sand","AD_sand","Watson_sand")
for(i in 1:M){
  sample = rnorm(n=n,mean=0,sd=1)
  mle = estimate.normal(sample)
  pval[i,1] = CvM.normal.pvalue(w=CvM.normal(x=sample))$P
  pval[i,2] = AD.normal.pvalue(a=AD.normal(x=sample))$P
  pval[i,3] = Watson.normal.pvalue(u=Watson.normal(x=sample))$P
  pval[i,4] = gof.sandwich(y=sample,Fdist=cdf.normal.user,thetahat=mle,Score=score.normal.user,m=100)$CvM$P
  pval[i,5] = gof.sandwich(y=sample,Fdist=cdf.normal.user,thetahat=mle,Score=score.normal.user,m=100)$AD$P
  pval[i,6] = gof.sandwich(y=sample,Fdist=cdf.normal.user,thetahat=mle,Score=score.normal.user,m=100)$Watson$P
}
par(mfrow=c(1,3))
plot(x=pval[,1],y=pval[,4],main="CvM",ylab="sandwich")
plot(x=pval[,2],y=pval[,5],main="AD")
plot(x=pval[,3],y=pval[,6],main="Watson")

par(mfrow=c(2,3))
hist(pval[,1],main="CvM")
hist(pval[,2],main="AD")
hist(pval[,3],main="Watson")
hist(pval[,4],ylab="sandwich")
hist(pval[,5])
hist(pval[,6])
