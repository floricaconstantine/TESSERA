library(testthat)
library(Matrix)
# Ensure sparseinv is loaded since E_step_Vhat depends on it
library(sparseinv)

# --- Setup Shared Test Data ---
set.seed(2026)
n <- 5
# Convert Q_test to a sparse dgCMatrix to mirror the spatial precision matrix
Q_test <- methods::as(Matrix::rsparsematrix(n, n, density = 0.5, symmetric = TRUE),
                      "dgCMatrix")
diag(Q_test) <- rowSums(abs(Q_test)) + 1
tau2_test <- 0.5
z_test <- rpois(n, 10)
X_test <- matrix(1, nrow = n, ncol = 1)
beta_test <- 0.2
N_test <- rep(100, n)

test_that("E_step_Vhat returns a valid subset matching the exact dense inverse diagonal",
          {
            Vhat_sparse <- E_step_Vhat(Q_test, tau2_test, z_test)
            
            # Ensure it returns a matrix-like object
            expect_true(inherits(Vhat_sparse, "Matrix") ||
                          is.matrix(Vhat_sparse))
            
            # Compute the exact dense inverse manually
            V_inv_manual <- (Q_test / tau2_test) + Matrix::Diagonal(n, 0.5 + z_test)
            V_exact_dense <- Matrix::solve(V_inv_manual)
            
            # The Takahashi subset MUST exactly match the diagonal of the true dense inverse
            expect_equal(as.numeric(Matrix::diag(Vhat_sparse)),
                         as.numeric(Matrix::diag(V_exact_dense)),
                         tolerance = 1e-8)
          })

test_that("E_step_Vhat handles inversion failure gracefully via fallbacks",
          {
            # Create a perfectly singular matrix to trigger failures
            Q_singular <- Matrix::Matrix(0, 3, 3, sparse = TRUE)
            z_singular <- rep(-0.5, 3) # Creates a 0 diagonal in Vinv
            
            # Suppress all the expected fallback warnings so they don't clutter the console,
            # and capture the final returned matrix.
            suppressWarnings({
              res <- E_step_Vhat(Q_singular, 1, z_singular)
            })
            
            # Verify that despite the internal failures, it successfully fell back to the
            # pseudoinverse and returned a valid matrix
            expect_true(inherits(res, "matrix") || inherits(res, "Matrix"))
            expect_equal(dim(res), c(3, 3))
          })

test_that("E_step_etahat computes the correct sparse linear system solution",
          {
            # Note: Vhat is no longer passed into E_step_etahat
            eta_hat <- E_step_etahat(Q_test, tau2_test, beta_test, X_test, z_test, N_test)
            
            # Manual reconstruction of the right-hand side
            term1 <- z_test + 0.5
            term1 <- term1 * log(term1 / N_test) - 0.5
            term2 <- (Q_test %*% (X_test %*% beta_test)) / tau2_test
            
            # Manual reconstruction of Vinv
            Vinv <- (Q_test / tau2_test) + Matrix::Diagonal(n, 0.5 + z_test)
            
            # The expectation should be the exact solution to Vinv * eta = (term1 + term2)
            expected_eta <- as.numeric(Matrix::solve(Vinv, term1 + term2))
            
            expect_equal(as.numeric(eta_hat), expected_eta, tolerance = 1e-8)
          })

test_that("E_step_thetahat handles log-linear transformation correctly with sparse diagonals",
          {
            # Mock a sparse covariance matrix subset
            Vhat <- Matrix::Diagonal(n, 0.1)
            eta_hat <- rep(0, n)
            
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
