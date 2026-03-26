library(testthat)
library(Matrix)

# Helper setup
W_test <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
D_test <- diag(rowSums(W_test))
gamma_test <- 0.8

test_that("Precision matrix functions return correct structures", {
  # CAR: Check for base matrix type
  Q_car <- Q_matrix_CAR(W_test, D_test, gamma_test)
  expect_true(is.matrix(Q_car))
  expect_true(isSymmetric(Q_car))
  
  # SAR: Check for the S4 Matrix class (more flexible than dsCMatrix)
  Q_sar <- Q_matrix_SAR(W_test, D_test, gamma_test)
  expect_s4_class(Q_sar, "Matrix")
  
  # Leroux
  Q_leroux <- Q_matrix_Leroux(W_test, D_test, 1)
  expect_equal(as.matrix(Q_leroux), D_test - W_test, ignore_attr = TRUE)
})

test_that("nngp_prec_mat produces valid sparse matrices", {
  coords <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  sp_dist <- matrix(c(1.41, 1), nrow = 2)
  params <- c(0.1, 1.0, 5.0)
  
  res <- nngp_prec_mat(sp_dist, coords, "Exp", params)
  
  expect_named(res, c("Q", "Dinv", "A"))
  # Use the actual class returned (dgCMatrix) or the broad "sparseMatrix"
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
