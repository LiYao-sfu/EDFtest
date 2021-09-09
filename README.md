# GOFEmpirical
This repository contains software for the calculation of goodness-of-fit
test statistics and their P-values.  The three statistics computed are the
Empirical Distribution function statistics called Cramér-von Mises, Anderson-Darling,
and Watson's statistic.  

The statistics and their P-values can be used to assess an assumed distribution. In the simplest situation
you have an iid sample from some distribution $F$ and want to test the hypothesis that $F$ is a member of 
some specific parametric family. The following families are available:

Normal(\mu,\sigma^2)
Gamma(shape = \alpha, scale = \beta)
Logistic(location = \mu, scale = \beta)
Laplace(location = \mu, scale = \beta)
Weibull(shape = \alpha, scale = \beta)
Extreme Value(location = \mu, scale = \beta)

The package also includes regression models in which a response $Y$ is related to predictors $X$. 
The model specifies the conditional distribution of $Y$ given $X$.  The package contains code
for situations where the conditional distribution is one of the list given above.  The 
following models are handled:

Linear regression with homosecdastic errors: Y_i has a N(X_i \beta, \sigma^2) distribution given X_i.
