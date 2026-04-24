library(testthat)
library(Matrix)

# --- Shared Setup for Lattice Wrapper Tests ---
set.seed(2026)
n_test_pts <- 100
nb_test_dist <- 0.15

# Generate two areas with clear spatial signal
area1 <- generate_data_one_area(
  n_points = n_test_pts,
  nb_dist = nb_test_dist,
  model_type = "Leroux",
  beta_true = c(1, -0.5),
  gamma_true = 0.5,
  tau2_true = 0.8
)
area2 <- generate_data_one_area(
  n_points = n_test_pts,
  nb_dist = nb_test_dist,
  model_type = "Leroux",
  beta_true = c(1, -0.5),
  gamma_true = 0.5,
  tau2_true = 0.8
)

mock_data <- structure(
  list(
    counts_list = list(
      A1 = matrix(area1$z, nrow = 1),
      A2 = matrix(area2$z, nrow = 1)
    ),
    X_list = list(A1 = area1$X, A2 = area2$X),
    W_list = list(A1 = area1$W, A2 = area2$W),
    D_list = list(A1 = area1$D, A2 = area2$D),
    library_size_list = list(A1 = area1$library_size, A2 = area2$library_size),
    coords_list = list(A1 = area1$coords, A2 = area2$coords),
    eig_L_list = list(A1 = area1$eig_list, A2 = area2$eig_list),
    eig_CS_list = NULL
  ),
  class = "TESSERAData"
)
rownames(mock_data$counts_list$A1) <- "Gene1"
rownames(mock_data$counts_list$A2) <- "Gene1"


# --- Tests ---

test_that("check_inputs_TESSERA validates object structure correctly", {
  expect_silent(check_inputs_TESSERA(mock_data))
  
  broken_data <- mock_data
  broken_data$W_list <- NULL
  expect_error(check_inputs_TESSERA(broken_data))
})

test_that("TESSERA_lattice fits multiple model types and handles eigs", {
  # Test 1: Leroux (Uses pre-computed eigs)
  res_l <- suppressMessages(suppressWarnings(
    TESSERA_lattice(
      mock_data,
      gene_name = "Gene1",
      model_type = "Leroux",
      em_iters = 3,
      verbose = FALSE
    )
  ))
  expect_s3_class(res_l, "TESSERAOutput")
  expect_length(res_l$beta_hat, 2)
  
  # Test 2: CAR (Triggers on-the-fly computation of CS eigenvalues)
  res_c <- suppressMessages(suppressWarnings(
    TESSERA_lattice(
      mock_data,
      gene_name = "Gene1",
      model_type = "CAR",
      em_iters = 3,
      verbose = FALSE
    )
  ))
  expect_named(res_c$gamma_hat, c("A1", "A2"))
})

test_that("TESSERA_lattice respects early stopping criteria", {
  res <- suppressMessages(suppressWarnings(
    TESSERA_lattice(
      mock_data,
      gene_name = "Gene1",
      model_type = "Leroux",
      em_iters = 20,
      em_min_iters = 2,
      em_tol = 1.0,
      em_stopping = "abs_beta_norm",
      verbose = FALSE
    )
  ))
  
  expect_true(res$run_settings$em_iters_actual < 20)
})

test_that("TESSERA_lattice produces valid trackers and Hessians", {
  res <- suppressMessages(suppressWarnings(
    TESSERA_lattice(
      mock_data,
      gene_name = "Gene1",
      em_iters = 5,
      verbose = FALSE
    )
  ))
  
  # Tracker dimensions
  expect_equal(ncol(res$beta_tracker), 6)
  
  # FIXED: beta_neghessian is an S4 Matrix object, not a base double vector/matrix
  expect_s4_class(res$beta_neghessian, "Matrix")
  
  # tau2_neghessian is a standard numeric vector, so "double" is correct here
  expect_type(res$tau2_neghessian, "double")
  expect_named(res$tau2_neghessian, c("A1", "A2"))
})
