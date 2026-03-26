library(testthat)
library(Matrix)

# --- Shared Setup for M-Step Tests ---
set.seed(2026)
n <- 10
Q_test <- diag(n)
Vhat_test <- diag(rep(0.1, n))
eta_test <- rnorm(n)
X_test <- matrix(1, nrow = n, ncol = 1)
beta_test <- 0.5

# --- M_step_tau2 ---
test_that("M_step_tau2 computes valid scaling parameters", {
  # 1. Standard calculation
  res <- M_step_tau2(Vhat_test, eta_test, Q_test, beta_test, X_test)
  
  expect_type(res, "double")
  expect_gt(res, 0)
  
  # 2. Test the "off the rails" failsafe
  Q_zero <- matrix(0, 2, 2)
  V_zero <- matrix(0, 2, 2)
  eta_zero <- c(0, 0)
  
  # FIX: Put the assignment INSIDE the expectation.
  # This ensures expect_warning captures the warning,
  # but allows res_small to be assigned the numeric return value.
  expect_warning(
    res_small <- M_step_tau2(V_zero, eta_zero, Q_zero, 0, matrix(0, 2, 1)),
    "Invalid tau^2",
    fixed = TRUE
  )
  
  # Now res_small will be the numeric value, not the warning object
  expect_equal(as.numeric(res_small), 100 * .Machine$double.eps)
})

# --- M_step_beta ---
test_that("M_step_beta aggregates across multiple areas correctly", {
  eta_list <- list(rep(1, 5), rep(1, 5))
  Q_list <- list(diag(5), diag(5))
  tau2_list <- list(1, 1)
  X_list <- list(matrix(1, 5, 1), matrix(1, 5, 1))
  
  res <- M_step_beta(eta_list, Q_list, tau2_list, X_list)
  expect_equal(as.numeric(res), 1.0)
})

# --- M_step_gamma_CAR ---
test_that("M_step_gamma_CAR finds a root within bounds", {
  eig_vals <- seq(0.1, 0.9, length.out = 10)
  W <- diag(10)
  D <- diag(10)
  
  res <- M_step_gamma_CAR(Vhat_test, eta_test, 1.0, beta_test, X_test, W, D, eig_vals)
  
  expect_named(res, c("gamma_hat", "grad_val"))
  expect_true(res$gamma_hat >= 0 && res$gamma_hat < 1)
})

# --- M_step_variogram ---
test_that("M_step_variogram returns expected parameter count", {
  skip_if_not_installed("gstat")
  skip_if_not_installed("sp")
  
  coords <- matrix(runif(20), ncol = 2)
  eta_hat <- rnorm(10)
  X <- matrix(1, 10, 1)
  beta <- 0
  
  # Suppress internal gstat "No convergence" warnings
  res <- suppressWarnings(M_step_variogram(eta_hat, beta, X, coords, cov_type = "Exp"))
  
  expect_length(res, 3)
  expect_true(all(res >= 0))
})

# --- M_step_BRISC ---
test_that("M_step_BRISC wrapper functions correctly", {
  skip_if_not_installed("BRISC")
  
  coords <- matrix(runif(20), ncol = 2)
  eta_hat <- rnorm(10)
  X <- matrix(1, 10, 1)
  beta <- 0
  
  # Suppress the failsafe warnings here too.
  # It fails because 10 random points rarely have a clear spatial structure,
  # but we want to verify the function still returns the failsafe numeric vector.
  res <- suppressWarnings(M_step_BRISC(eta_hat, beta, X, coords, k = 5))
  
  expect_type(res, "double")
  expect_true(length(res) >= 3)
})
