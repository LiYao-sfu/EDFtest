# Test functions: CvM.normal.covmat, CvM.logistic.covmat, CvM.weibull.covmat
# Function missing: Gamma

#Used by p-values functions, when p value uniform

n=10
CvM.normal.covmat(n) #typo in the function
CvM.weibull.covmat(n)
CvM.logistic.covmat(n)

#Ask Richard: How to test these three functions

# Test functions: AD.normal.covmat, AD.logistic.covmat, AD.weibull.covmat

AD.normal.covmat(n)
AD.weibull.covmat(n)
AD.logistic.covmat(n)


# Test functions: CvM.logistic.eigen, CvM.weibull.eigen
# Function missing: normal, gamma

CvM.logistic.eigen(n)
CvM.weibull.eigen(n)
AD.logistic.eigen(n)
AD.weibull.eigen(n)


