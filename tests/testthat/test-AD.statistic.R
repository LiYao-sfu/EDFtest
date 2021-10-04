context("test-AD.statistic.R")

test_that("AD for normal sample", {
  x = c(0.25024690, -0.33712454, -0.11335370, -0.09888291, 0.26408682,
        0.13898369, -0.24226950, 0.05903138, -0.17727187, 0.79468027)
  asq = AD.normal(x)
  asq_par = AD.normal(x,c(0,1))

  expect_equal(asq,0.380504793)
  expect_equal(asq_par,1.78788894)
})

test_that("AD for gamma sample", {
  x = c(0.5047757, 0.1538300, 0.5704100, 0.3013008, 1.2775724,
        1.0468233, 0.6525627, 0.4376768, 2.4700737, 1.0944885)
  asq = AD.gamma(x)
  asq_par = AD.gamma(x,c(1,1))

  expect_equal(asq,0.201390974)
  expect_equal(asq_par,0.46290117)
})

test_that("AD for logistic sample", {
  x = c(-1.1263960,  0.9562103, -3.3860294,  0.1980448,  0.7667096,
        -0.8461510, -0.4524666,  1.0070690,  3.2450939,  1.1559508)
  asq = AD.logistic(x)
  asq_par = AD.logistic(x,c(0,1))

  expect_equal(asq,0.27184557)
  expect_equal(asq_par,0.35324973)
})

test_that("AD for laplace sample", {
  x = c(-0.23539279,  0.16009027,  2.84634962,  0.35710312, -0.40466195,
        -0.41113889,  2.16169132, -0.27151351,  0.13770907,  0.02330074)
  asq = AD.laplace(x)
  asq_par = AD.laplace(x,c(0,1))

  expect_equal(asq,0.91031712)
  expect_equal(asq_par,1.0260448)
})

test_that("AD for weibull sample", {
  x = c(0.36218715, 0.16506700, 0.16757965, 0.93681048, 1.87396510,
        0.44718470, 1.24767735, 0.07435952, 1.86023456, 0.03682825)
  asq = AD.weibull(x)
  asq_par = AD.weibull(x,c(1,1))

  expect_equal(asq,0.30993253)
  expect_equal(asq_par,0.8129256)
})

test_that("AD for exponential sample", {
  x = c(13.0581121,  0.8301048,  0.5207504,  1.0923122,  0.7086793,
        0.1271974,  3.9326089,  0.0510448,  4.3839846,  3.4396530)
  asq = AD.exp(x)
  asq_par = AD.exp(x,2)

  expect_equal(asq,0.93011782)
  expect_equal(asq_par,0.8513746)
})
