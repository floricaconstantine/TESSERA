library(testthat)
library(Matrix)

# Helper setup
W_test <- Matrix::Matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3, sparse = TRUE)
D_test <- Matrix::Diagonal(3, rowSums(as.matrix(W_test)))
gamma_test <- 0.8

# --- Precomputations for new signatures ---
# SAR: W D^{-1} W
D_inv_test <- Matrix::Diagonal(n = nrow(D_test), x = 1 / Matrix::diag(D_test))
precomp_sar_test <- W_test %*% D_inv_test %*% W_test

# Leroux: list(D_minus_W, id_mat)
precomp_leroux_test <- list(D_minus_W = D_test - W_test,
                            id_mat = Matrix::Diagonal(n = nrow(W_test), 1))
# ----------------------------------------

test_that("Precision matrix functions return correct structures", {
  # CAR: Check for base matrix type (precomp not needed)
  Q_car <- Q_matrix_CAR(W_test, D_test, gamma_test)
  # Q_car could be sparse Matrix depending on inputs, so we check inheritence
  expect_true(inherits(Q_car, "Matrix") || is.matrix(Q_car))
  expect_true(isSymmetric(as.matrix(Q_car)))
  
  # SAR: Check for the S4 Matrix class and correct mathematical expansion
  Q_sar <- Q_matrix_SAR(W_test, D_test, gamma_test, precomp = precomp_sar_test)
  expect_s4_class(Q_sar, "Matrix")
  
  # Mathematical verification of SAR expansion
  Q_sar_expected <- D_test - (2 * gamma_test * W_test) + ((gamma_test^2) * precomp_sar_test)
  expect_equal(as.matrix(Q_sar), as.matrix(Q_sar_expected), ignore_attr = TRUE)
  
  # Leroux: Check with gamma = 1
  Q_leroux_gamma1 <- Q_matrix_Leroux(W_test, D_test, 1, precomp = precomp_leroux_test)
  expect_equal(as.matrix(Q_leroux_gamma1),
               as.matrix(D_test - W_test),
               ignore_attr = TRUE)
  
  # Leroux: Check with gamma = 0 to ensure identity matrix precomp triggers correctly
  Q_leroux_gamma0 <- Q_matrix_Leroux(W_test, D_test, 0, precomp = precomp_leroux_test)
  expect_equal(
    as.matrix(Q_leroux_gamma0),
    as.matrix(precomp_leroux_test$id_mat),
    ignore_attr = TRUE
  )
})

test_that("nngp_prec_mat produces valid sparse matrices with precomputed distances",
          {
            # sp_dist has 1 column.
            # Row 1 is the distance (1.41), Row 2 is the neighbor index (1).
            sp_dist <- matrix(c(1.41, 1), nrow = 2)
            
            # Mock the new precomputed nb_dist list
            # Since there is only 1 neighbor, the distance matrix between neighbors is 1x1 zero.
            nb_dist <- list(matrix(0.0, 1, 1))
            
            params <- c(0.1, 1.0, 5.0)
            
            # Signature updated to pass nb_dist instead of coords
            res <- nngp_prec_mat(sp_dist, nb_dist, "Exp", params)
            
            expect_named(res, c("Q", "Dinv", "A"))
            expect_s4_class(res$Q, "sparseMatrix")
            expect_length(res$Dinv, 2)
          })

test_that("expected_loglike catches invalid model types", {
  expect_error(
    expected_loglike(
      diag(2),
      c(1, 1),
      diag(2),
      0.5,
      1,
      c(1, 1),
      matrix(1, 2, 1),
      diag(2),
      diag(2),
      c(1, 1),
      "INVALID"
    ),
    "Invalid model_type"
  )
})

test_that("neg_hessian_tau2 computes correctly using sparse extraction", {
  n <- 3
  X <- matrix(1, n, 1)
  beta <- 0.5
  eta <- c(0.1, -0.2, 0.3)
  tau2 <- 1.5
  Q <- methods::as(Q_matrix_CAR(W_test, D_test, gamma_test), "dgCMatrix")
  
  # Mock a sparse Vhat subset
  Vhat <- methods::as(Matrix::solve(Q / tau2 + Matrix::Diagonal(n, 1)), "dgCMatrix")
  
  res <- neg_hessian_tau2(Vhat, eta, Q, tau2, beta, X)
  
  # Manual dense computation to verify equivalence
  vector_term <- eta - (X %*% beta)
  term1 <- as.numeric(crossprod(Q %*% vector_term, vector_term))
  term2 <- sum(as.matrix(Q) * as.matrix(Vhat)) # Dense trace trick
  expected <- (term1 + term2) / (tau2^3) - (0.5 * nrow(X)) / (tau2^2)
  
  expect_equal(res, expected)
})

test_that("neg_hessian_gamma_SAR computes correctly with WZ sparse extraction",
          {
            n <- 3
            X <- matrix(1, n, 1)
            beta <- 0.5
            eta <- c(0.1, -0.2, 0.3)
            tau2 <- 1.5
            Q <- methods::as(Q_matrix_SAR(W_test, D_test, gamma_test), "dgCMatrix")
            Vhat <- methods::as(Matrix::solve(Q / tau2 + Matrix::Diagonal(n, 1)), "dgCMatrix")
            
            # Symmetric eigenvalues for SAR
            D_inv_sqrt <- Matrix::Diagonal(n, 1 / sqrt(Matrix::diag(D_test)))
            S <- D_inv_sqrt %*% W_test %*% D_inv_sqrt
            eig_vals <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
            
            res <- neg_hessian_gamma_SAR(Vhat, eta, gamma_test, tau2, beta, X, W_test, D_test, eig_vals)
            
            # Simply checking that the sparse index extraction runs and returns a valid scalar
            expect_type(res, "double")
            expect_length(res, 1)
            expect_true(is.finite(res))
          })

test_that("expected_loglike computes correct value with sparse Vhat extraction",
          {
            n <- 3
            X <- matrix(1, n, 1)
            beta <- 0.5
            eta <- c(0.1, -0.2, 0.3)
            tau2 <- 1.5
            Q <- methods::as(Q_matrix_CAR(W_test, D_test, gamma_test), "dgCMatrix")
            Vhat <- methods::as(Matrix::solve(Q / tau2 + Matrix::Diagonal(n, 1)), "dgCMatrix")
            
            D_inv_sqrt <- Matrix::Diagonal(n, 1 / sqrt(Matrix::diag(D_test)))
            S <- D_inv_sqrt %*% W_test %*% D_inv_sqrt
            eig_vals <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
            
            res <- expected_loglike(Vhat,
                                    eta,
                                    Q,
                                    gamma_test,
                                    tau2,
                                    beta,
                                    X,
                                    W_test,
                                    D_test,
                                    eig_vals,
                                    "CAR")
            
            expect_type(res, "double")
            expect_length(res, 1)
            expect_true(is.finite(res))
          })
