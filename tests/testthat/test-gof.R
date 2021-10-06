context("test-gof.R")

test_that("gof imhof method for normal sample", {
  x = c(0.25024690, -0.33712454, -0.11335370, -0.09888291, 0.26408682,
        0.13898369, -0.24226950, 0.05903138, -0.17727187, 0.79468027)
  result = gof.normal(x)

  expect_equal(result$Wsq,0.051555557)
  expect_equal(result$Wsq.pvalue,0.49021528)
  expect_equal(result$Asq,0.380504793)
  expect_equal(result$Asq.pvalue,0.39529586)
  expect_equal(result$Usq,0.044353005)
  expect_equal(result$Usq.pvalue,0.56491446)

  expect_output(str(result), "List of 6")
})

test_that("gof bootstrap method for normal sample", {
  x = c(0.25024690, -0.33712454, -0.11335370, -0.09888291, 0.26408682,
        0.13898369, -0.24226950, 0.05903138, -0.17727187, 0.79468027)
  result = gof.normal.bootstrap(x,M=100)

  expect_equal(result$Wsq,0.051555557)
  expect_equal(result$Asq,0.380504793)
  expect_equal(result$Usq,0.044353005)

  expect_output(str(result), "List of 6")
})

test_that("gof imhof method for gamma sample", {
  x = c(0.5047757, 0.1538300, 0.5704100, 0.3013008, 1.2775724,
        1.0468233, 0.6525627, 0.4376768, 2.4700737, 1.0944885)
  result = gof.gamma(x)

  expect_equal(result$Wsq,0.029695982)
  expect_equal(result$Wsq.pvalue,0.49021528) #start from here!
  expect_equal(result$Asq,0.380504793)
  expect_equal(result$Asq.pvalue,0.39529586)
  expect_equal(result$Usq,0.044353005)
  expect_equal(result$Usq.pvalue,0.56491446)

  expect_output(str(result), "List of 6")
})

test_that("gof bootstrap method for gamma sample", {
  x = c(0.5047757, 0.1538300, 0.5704100, 0.3013008, 1.2775724,
        1.0468233, 0.6525627, 0.4376768, 2.4700737, 1.0944885)
  result = gof.gamma.bootstrap(x,M=100)

  expect_equal(result$Wsq,0.029695982)
  expect_equal(result$Asq,0.380504793)
  expect_equal(result$Usq,0.044353005)

  expect_output(str(result), "List of 6")
})
