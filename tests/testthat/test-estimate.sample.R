context("test-estimate.sample.R")

test_that("mle estimates for normal distribution", {
  x = rnorm(n=100,mean=0,sd=1)
  par = estimate.normal(x)

  expect_equal(length(par),2)
  #mean
  expect_lt(par[1],qnorm(0.9999)/sqrt(100))
  expect_gt(par[1],qnorm(0.0001)/sqrt(100))
  #sd
  expect_lt(par[2],sqrt(qchisq(0.9999,df=99)/99))
  expect_gt(par[2],sqrt(qchisq(0.0001,df=99)/99))
})

test_that("mle estimates for gamma distribution", {
  set.seed(100)
  x = rgamma(n=100,shape=1,scale=2)
  par.scale = estimate.gamma(x)
  par.rate = estimate.gamma(x,use.rate = TRUE)

  expect_equal(length(par.scale),2)

  expect_lt(par.scale[1],)
  expect_gt(par.scale[1],0.)

  expect_lt(par.scale[2],2.5)
  expect_gt(par.scale[2],1.5)

  expect_equal(par.scale[2],1/par.rate[2])
})

test_that("mle estimates for logistic distribution", {
  x = rlogis(n=100,location=0,scale=1)
  par = estimate.logistic(x,verbose = TRUE)

  expect_lt(par[1],0.5)
  expect_gt(par[1],-0.5)

  expect_lt(par[2],1.5)
  expect_gt(par[2],0.5)
})

test_that("mle estimates for laplace distribution", {
  x1 = rmutil::rlaplace(n=100,m=0,s=1)
  x2 = L1pack::rlaplace(n=100,location=0,scale=1)
  par.sd = estimate.laplace(x1)
  par.scale = estimate.laplace(x2,use.sd=FALSE)

  expect_equal(par.sd,par.scale)

})

test_that("mle estimates for weibull distribution", {
  x = rweibull(n=100,shape=1,scale=1)
  par = estimate.weibull(x)

  expect_lt(par[1],1.5)
  expect_gt(par[1],0.5)

  expect_lt(par[2],1.5)
  expect_gt(par[2],0.5)
})

test_that("mle estimates for exponential distribution", {
  x = rexp(n=100,rate=1/2)
  par.rate = estimate.exp(x,use.rate=TRUE)
  par.scale = estimate.exp(x,use.rate=FALSE)

  expect_equal(par.rate,1/par.scale)

  expect_lt(par.scale,2.5)
  expect_gt(par.scale,1.5)

})

