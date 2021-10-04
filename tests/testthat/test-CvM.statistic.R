context("test-CvM.statistic.R")

test_that("CvM statistics for normal sample", {
  wsq <- rep(0,100)
  for(i in 1:100){
    x = rnorm(n=100,mean=0,sd=1)
    wsq[i] = CvM.normal(x)
    CvM.normal(x,c(0,1))
  }
  avg = mean(wsq)
  expect_lt(avg,0.08)
  expect_gt(avg,0)
})

test_that("CvM statistics for gamma sample", {
  x = rgamma(n=100,shape=1,scale=2)
  CvM.gamma(x)
  CvM.gamma(x,c(0,1))
})
