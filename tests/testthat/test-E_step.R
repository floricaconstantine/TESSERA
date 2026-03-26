library(testthat)
library(Matrix)

# --- Setup Shared Test Data ---
set.seed(2026)
n <- 5
Q_test <- as.matrix(Matrix::rsparsematrix(n, n, density = 0.5, symmetric = TRUE))
diag(Q_test) <- rowSums(abs(Q_test)) + 1
tau2_test <- 0.5
z_test <- rpois(n, 10)
X_test <- matrix(1, nrow = n, ncol = 1)
beta_test <- 0.2
N_test <- rep(100, n)

test_that("E_step_Vhat and E_step_Vhat_PLU are consistent", {
  Vhat <- E_step_Vhat(Q_test, tau2_test, z_test)
  decomp <- E_step_Vhat_PLU(Q_test, tau2_test, z_test)
  
  # FIX 1: Check names regardless of their order in the list
  expect_setequal(names(decomp), c("P", "L", "U", "Linv", "Uinv"))
  
  # Verify PLU reconstruction
  V_inv_manual <- (Q_test / tau2_test) + diag(0.5 + z_test)
  reconstructed <- as.matrix(decomp$P %*% decomp$L %*% decomp$U)
  expect_equal(reconstructed, as.matrix(V_inv_manual), ignore_attr = TRUE)
})

test_that("E_step_Vhat handles inversion failure gracefully via pinv", {
  # Create a singular matrix
  Q_singular <- matrix(0, 3, 3)
  z_singular <- rep(-0.5, 3)
  
  # FIX 2: Since your function uses pracma::pinv() as a fallback,
  # it likely generates a warning and returns a matrix rather than an error.
  expect_warning({
    res <- E_step_Vhat(Q_singular, 1, z_singular)
  }, "INVERSION OF Inv\\(V_hat\\) FAILED")
  
  # Verify that it still returned a matrix (the pseudoinverse)
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(3, 3))
})

test_that("E_step_etahat results match between standard and PLU methods", {
  Vhat <- E_step_Vhat(Q_test, tau2_test, z_test)
  decomp <- E_step_Vhat_PLU(Q_test, tau2_test, z_test)
  
  # Standard calculation
  eta_std <- E_step_etahat(Vhat, Q_test, tau2_test, beta_test, X_test, z_test, N_test)
  
  # PLU calculation
  eta_plu <- E_step_etahat_PLU(decomp, Q_test, tau2_test, beta_test, X_test, z_test, N_test)
  
  expect_equal(as.numeric(eta_std), as.numeric(eta_plu), tolerance = 1e-8)
})

test_that("E_step_thetahat handles log-linear transformation correctly",
          {
            Vhat <- diag(rep(0.1, n)) # Mock covariance
            eta_hat <- rep(0, n)       # Mock mean
            
            # theta = exp(eta + 0.5 * diag(V))
            expected_theta <- exp(0 + 0.5 * 0.1)
            res <- E_step_thetahat(Vhat, eta_hat)
            
            expect_equal(as.numeric(res), rep(expected_theta, n))
          })

test_that("E_step_predict follows the Poisson-Gaussian back-transformation",
          {
            theta_hat <- c(0.5, 1.0)
            N <- c(100, 100)
            # formula: theta * N - 0.5
            expected <- (theta_hat * N) - 0.5
            
            res <- E_step_predict(theta_hat, N)
            expect_equal(res, expected)
          })
