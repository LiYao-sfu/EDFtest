context("test-CvM.pvalue.R")

test_that("CvM P-value for normal sample", {
  pval <- rep(0,100)
  for(i in 1:100){
    x = rnorm(n=100,mean=0,sd=1)
    wsq = CvM.normal(x)
    pval[i] = CvM.normal.pvalue(wsq)$P
  }
  avg = mean(pval)
  expect_lt(avg,0.4)
  expect_gt(avg,0.6)
})
