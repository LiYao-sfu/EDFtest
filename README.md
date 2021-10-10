# EDFtest
This repository contains software for the calculation of goodness-of-fit
test statistics and their P-values.  The three statistics computed are the
Empirical Distribution function statistics called Cramér-von Mises, Anderson-Darling,
and Watson statistic.  

The statistics and their P-values can be used to assess an assumed distribution. In the simplest situation
you have an i.i.d. sample from some distribution F and want to test the hypothesis that the sample is drawn from
a distribution F which belongs to a specified parametric family of distributions against the alternative that 
F is not equal to any member of that parametric family. The following families are available:
Uniform(min,max)
Normal(location,scale),
Gamma(shape,scale),
Logistic(location,scale),
Laplace(location,scale),
Weibull(shape,scale), and
Exponential(scale).

Users can add their own distributions by providing two functions

* `Fdist(x,thetahat, ...)` which takes parameter estimates and a data set `x` and computes the probability integral transform for each element of `x`. `Fdist` must return a vector of n probabilities

* `Score(x,thetahat, ...)` which takes parameter estimates and a data set `x` and computes, for each entry in `x`, the component of the score function due to observation `x`. These must be returned in an $n$ by $p$ matrix with 1 row for each observation and 1 column for each parameter.

The user is also expected to supply the value `thetahat` of the maximum likelihood estimate.

The package also includes regression models in which a response Y is related to predictors X. 
The model specifies the conditional distribution of Y given X.  The package contains code
for situations where the conditional distribution is one of the list given above.  The 
following models are handled:

Linear regression with homosecdastic errors: <img src="https://render.githubusercontent.com/render/math?math=Y_i"> has a N(<img src="https://render.githubusercontent.com/render/math?math=X_i \beta, \sigma^2">) distribution given <img src="https://render.githubusercontent.com/render/math?math=X_i">.

Authors:

-   [Li Yao](https://github.com/LiYao-sfu),
    <yaoliy@sfu.ca> (Maintainer)
-   [Richard Lockhart](http://www.sfu.ca/~lockhart/),
    <lockhart@sfu.ca>

Papers:

-   [Paper Title 1](https:) *Journal a*



## Installation
There are several ways you can install GitHub packages into R. For example,
You can install our package by using `devtools`. You need to install `devtools` package first if you have not.


Step 1: Install the `devtools` package
```R
install.packages("devtools")
```

Step 2: Install our `EDFtest` package and load it
```R
library(devtools)
install_github("LiYao-sfu/EDFtest")
library("EDFtest")
```

## Troubleshooting
This package is still under development. EDF test for regression models and discrete discrete distributions 
will be available for the next minor release.

If you encounter a clear bug, please create an issue on github. For questions and other discussion, please 
contact Li Yao by his email.
