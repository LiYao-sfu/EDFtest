context("test-Watson.statistic.R")

test_that("Watson for uniform sample", {
  x = c(0.3103205, -0.1386720, -0.5988813,  0.9483934, -0.5213117,
        0.7555062, -0.6821612,  0.7957394, -0.7387161,  0.1910647)
  usq = Watson.uniform(x)

  expect_equal(usq,0.1168008)
})

test_that("Watson for normal sample", {
  x = c(0.25024690, -0.33712454, -0.11335370, -0.09888291, 0.26408682,
        0.13898369, -0.24226950, 0.05903138, -0.17727187, 0.79468027)
  usq = Watson.normal(x)

  expect_equal(usq,0.044353005)
})

test_that("Watson for gamma sample", {
  x = c(0.5047757, 0.1538300, 0.5704100, 0.3013008, 1.2775724,
        1.0468233, 0.6525627, 0.4376768, 2.4700737, 1.0944885)
  usq = Watson.gamma(x)

  expect_equal(usq,0.027784624)
})

test_that("Watson for logistic sample", {
  x = c(-1.1263960,  0.9562103, -3.3860294,  0.1980448,  0.7667096,
        -0.8461510, -0.4524666,  1.0070690,  3.2450939,  1.1559508)
  usq = Watson.logistic(x)

  expect_equal(usq,0.039618567)
})

test_that("Watson for laplace sample", {
  x = c(-0.23539279,  0.16009027,  2.84634962,  0.35710312, -0.40466195,
        -0.41113889,  2.16169132, -0.27151351,  0.13770907,  0.02330074)
  usq = Watson.laplace(x)

  expect_equal(usq,0.131789187)
})

test_that("Watson for weibull sample", {
  x = c(0.36218715, 0.16506700, 0.16757965, 0.93681048, 1.87396510,
        0.44718470, 1.24767735, 0.07435952, 1.86023456, 0.03682825)
  usq = Watson.weibull(x)

  expect_equal(usq,0.04571839)
})

test_that("Watson for exponential sample", {
  x = c(13.0581121,  0.8301048,  0.5207504,  1.0923122,  0.7086793,
        0.1271974,  3.9326089,  0.0510448,  4.3839846,  3.4396530)
  usq = Watson.exp(x)

  expect_equal(usq,0.088244806)
})

