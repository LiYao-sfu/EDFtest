# GOFEmpirical
This repository contains software for the calculation of goodness-of-fit
test statistics and their P-values.  The three statistics computed are the
Empirical Distribution function statistics called Cramér-von Mises, Anderson-Darling,
and Watson's statistic.  

The statistics and their P-values can be used to assess an assumed distribution. In the simplest situation
you have an iid sample from some distribution <img src="https://render.githubusercontent.com/render/math?math=F"> and want to test the hypothesis that <img src="https://render.githubusercontent.com/render/math?math=F"> is a member of 
some specific parametric family. The following families are available:

<img src="https://render.githubusercontent.com/render/math?math=Normal(\mu,\sigma^2)">,
<img src="https://render.githubusercontent.com/render/math?math=Gamma(shape = \alpha, scale = \beta)">,
<img src="https://render.githubusercontent.com/render/math?math=Logistic(location = \mu, scale = \beta)">,
<img src="https://render.githubusercontent.com/render/math?math=Laplace(location = \mu, scale = \beta)">,
<img src="https://render.githubusercontent.com/render/math?math=Weibull(shape = \alpha, scale = \beta)">, and
<img src="https://render.githubusercontent.com/render/math?math=Extreme Value(location = \mu, scale = \beta)">

The package also includes regression models in which a response <img src="https://render.githubusercontent.com/render/math?math=Y"> is related to predictors <img src="https://render.githubusercontent.com/render/math?math=X">. 
The model specifies the conditional distribution of <img src="https://render.githubusercontent.com/render/math?math=Y"> given <img src="https://render.githubusercontent.com/render/math?math=X">.  The package contains code
for situations where the conditional distribution is one of the list given above.  The 
following models are handled:

Linear regression with homosecdastic errors: Y_i has a N(X_i \beta, \sigma^2) distribution given X_i.


# Installation
There are several ways you can install GitHub packages into R. For example,
You can install our package by using `devtools`. You need to install `devtools` package in R if you have not.


Step 1: Install the `devtools` package
```R
install.packages("devtools")
```

Step 2: Install our `GOFEmpirical` package and load it
```R
library(devtools)
install_github("EigenEye/GOF1.3/GOFEmpirical")
library("GOFEmpirical")
```

# Troubleshooting
