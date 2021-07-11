# Test functions: AD, CvM

library(GOFEmpirical)


#######################################################################

# Test for AD.R

m=10000
out <- matrix(0,m,1)
for(i in 1:m){
  z=runif(100,0,1)
  out[i,]=AD(z)
}
mean(out) #could be greater than 1 or close to 1

for(i in 1:m){
  z=runif(10,0,10)
  out[i,]=AD(z)
}
mean(out) #crush, since z is not bounded by 0 and 1

for(i in 1:m){
  z=runif(10,-1,1)
  out[i,]=AD(z)
}
mean(out)

#Comment: need to check whether input is valid or not

for(i in 1:m){
  z=rbeta(100,0.2,0.7)
  out[i,]=AD(z)
}
mean(out) #much greater than 1


#######################################################################

# Test for CvM.R

m=10000
out <- matrix(0,m,1)
for(i in 1:m){
  z=runif(100,0,1)
  out[i,]=CvM(z)
}
mean(out) #could be greater than 1/6 or close to 1/6

for(i in 1:m){
  z=runif(100,0,10)
  out[i,]=CvM(z)
}
mean(out) #not crush, but insensible result, large number, since z is not bounded by 0 and 1

for(i in 1:m){
  z=runif(100,-1,1)
  out[i,]=CvM(z)
}
mean(out) #not crush, but insensible result, large number, since z is not bounded by 0 and 1

for(i in 1:m){
  z=rbeta(100,0.2,0.7)
  out[i,]=AD(z)
}
mean(out) #much greater than 1, numbers close to 0 will destroy the statistics
