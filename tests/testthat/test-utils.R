library(testthat)

# --- moran_I_nb ---
test_that("moran_I_nb computes correct spatial autocorrelation", {
  # Simple case: 4 points in a line, clear spatial pattern
  # 1-2-3-4 where 1,2 are high and 3,4 are low
  y <- c(10, 10, 1, 1)
  W <- matrix(c(0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0),
              nrow = 4,
              byrow = TRUE)
  
  result <- moran_I_nb(y, W)
  
  expect_type(result, "double")
  expect_length(result, 1)
  # With this pattern, Moran's I should be positive
  expect_gt(result, 0)
})

# --- poisson_loglike ---
test_that("poisson_loglike handles valid and invalid numeric inputs", {
  z <- c(1, 2, 0)
  theta <- c(1, 2, 1)
  
  # Manual calculation: 1*log(1)-1 + 2*log(2)-2 + 0*log(1)-1
  # = (0-1) + (1.386-2) + (0-1) = -3.6137
  expected <- sum(pmax(0, z) * log(pmax(0, theta)) - pmax(0, theta))
  
  expect_equal(poisson_loglike(z, theta), expected)
  
  # Test the zeroing out of Inf/NaN as noted in the source
  expect_equal(poisson_loglike(c(0, 1), c(0, 1)), -1) # 0*log(0) is NaN, should be 0
})

# --- variables_from_list ---
test_that("variables_from_list instantiates variables in the target environment",
          {
            test_env <- new.env()
            my_list <- list(alpha = 10,
                            beta = "test",
                            gamma = c(1, 2, 3))
            
            variables_from_list(my_list, target_environ = test_env)
            
            expect_true(exists("alpha", envir = test_env))
            expect_equal(test_env$alpha, 10)
            expect_equal(test_env$beta, "test")
            expect_equal(test_env$gamma, c(1, 2, 3))
          })

# --- sparseDist & sparseDist_LT ---
test_that("sparseDist functions return correct dimensions and distances", {
  # Coordinates for a 3-4-5 triangle
  # A(0,0), B(3,0), C(0,4)
  coords <- matrix(c(0, 0, 3, 0, 0, 4), nrow = 3, byrow = TRUE)
  k <- 1
  
  # sparseDist (Upper Triangular)
  res_ut <- sparseDist(coords, k)
  expect_equal(dim(res_ut), c(2 * k, 2))
  
  # Use as.numeric() to drop the "d1" name attribute
  expect_equal(as.numeric(res_ut[1, 1]), 3)
  
  # sparseDist_LT (Lower Triangular)
  res_lt <- sparseDist_LT(coords, k)
  expect_equal(dim(res_lt), c(2 * k, 2))
  
  # Use as.numeric() to drop the "d1" name attribute
  expect_equal(as.numeric(res_lt[1, 2]), 4)
})
