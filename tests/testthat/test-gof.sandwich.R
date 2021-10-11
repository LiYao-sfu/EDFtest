context("test-gof.sandwich.R")

test_that("Statistics and P-value by sandwich estimator", {
  sample = c( 0.81878866, -0.03791044, -0.83140426, -0.51850492,  0.45387155,
             -1.28909386,  0.18681088, -2.11597297,  0.78370780,  0.37041184,
              2.00560921, -0.96623307, -1.16253597,  0.74950590, -0.36028086,
              0.70051008,  1.01708059, -0.66269175, -0.46364974,  2.19508727)
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
  output = gof.sandwich(y=sample,Fdist=cdf.normal.user,thetahat=mle,Score=score.normal.user,m=100)
  output

  expect_equal(output$CvM$W2,0.03635126)
  expect_equal(output$CvM$P,0.2083474)
  expect_equal(output$AD$A2,0.268198)
  expect_equal(output$AD$P,0.43739209)
  expect_equal(output$Watson$U2,0.43739209)
  expect_equal(output$CvM$P,0.43739209)
})
