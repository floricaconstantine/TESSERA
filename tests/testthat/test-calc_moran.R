library(testthat)

# --- Shared Setup ---
set.seed(2026)

# --- distanceCalculate ---

test_that("distanceCalculate computes inverse Euclidean distance", {
  # Distance between (0,0) and (3,4) is 5. Output should be 1/5 = 0.2
  expect_equal(distanceCalculate(0, 0, 3, 4), 0.2)
  
  # Distance to self is 0. Function should return 0 (not Inf) per code logic
  expect_equal(distanceCalculate(1, 1, 1, 1), 0)
})

# --- normalize ---

test_that("normalize centers a vector around its mean", {
  x <- c(1, 2, 3) # Mean is 2
  res <- normalize(x)
  
  expect_equal(res, c(-1, 0, 1))
  expect_equal(mean(res), 0)
})

# --- calc_moran ---

test_that("calc_moran returns valid Moran's I statistics", {
  n <- 10
  x <- rnorm(n)
  c1 <- runif(n)
  c2 <- runif(n)
  
  res <- calc_moran(x, c1, c2)
  
  # Result should be a numeric vector of length 3 (I, EI, SD)
  expect_length(res, 3)
  
  # Check Expected Value (EI) logic: EI = -1/(N-1)
  expected_ei <- -1 / (n - 1)
  expect_equal(res[2], expected_ei)
  
  # Standard deviation (SD) should be non-negative
  expect_true(res[3] >= 0)
})

test_that("calc_moran handles constant vectors or small samples", {
  # Constant vector results in zero denominator and NaN Moran's I
  # This verifies the function doesn't crash on edge cases
  res_const <- suppressWarnings(calc_moran(rep(5, 5), runif(5), runif(5)))
  expect_true(is.nan(res_const[1]))
})
