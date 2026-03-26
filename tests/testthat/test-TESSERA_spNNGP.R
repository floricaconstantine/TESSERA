library(testthat)
library(Matrix)

# --- Shared Setup for spNNGP Wrapper Tests ---
set.seed(2026)
n_test_pts <- 50 
cov_params_true <- c(0.1, 1.0, 0.2) 

# Generate mock spNNGP data
area1 <- generate_data_one_area_spNNGP(
  n_points = n_test_pts, nb_dist = 0.5, cov_type = "Exp", 
  cov_params = cov_params_true, nngp_k = 5, beta_true = c(1, -0.5)
)

mock_sp_data <- structure(
  list(
    counts_list = list(S1 = matrix(area1$z, nrow = 1)),
    X_list = list(S1 = area1$X),
    coords_list = list(S1 = area1$coords),
    library_size_list = list(S1 = area1$library_size)
  ),
  class = "TESSERAData"
)
rownames(mock_sp_data$counts_list$S1) <- "Gene1"


# --- Tests ---

test_that("checkInputsTESSERAspNNGP validates required spNNGP fields", {
  expect_silent(checkInputsTESSERAspNNGP(mock_sp_data))
  
  broken_data <- mock_sp_data
  broken_data$coords_list <- NULL
  expect_error(checkInputsTESSERAspNNGP(broken_data))
})

test_that("TESSERA_spNNGP fits using variogram method", {
  skip_if_not_installed("gstat")
  
  res <- suppressMessages(suppressWarnings(
    TESSERA_spNNGP(
      mock_sp_data, gene_name = "Gene1", cov_type = "Exp", 
      nngp_k = 5, em_iters = 3, cov_fit_method = "variogram", 
      verbose = FALSE
    )
  ))
  
  expect_s3_class(res, "TESSERAOutput")
  
  # Handle the dimension-dropping quirk: 
  # If it's a single area, it might be a vector (length 4) or matrix (1x4)
  param_hat <- res$cov_param_hat
  if (is.null(dim(param_hat))) {
    expect_length(param_hat, 4)
  } else {
    expect_equal(dim(param_hat), c(1, 4))
  }
  
  expect_named(res$beta_hat, colnames(mock_sp_data$X_list$S1))
})

test_that("TESSERA_spNNGP fits using BRISC method", {
  skip_if_not_installed("BRISC")
  
  res <- suppressMessages(suppressWarnings(
    TESSERA_spNNGP(
      mock_sp_data, gene_name = "Gene1", cov_init = "BRISC",
      cov_fit_method = "BRISC", em_iters = 2, verbose = FALSE
    )
  ))
  
  expect_true(res$run_settings$cov_fit_method == "BRISC")
  expect_s4_class(res$beta_neghessian, "Matrix")
})

test_that("TESSERA_spNNGP handles early stopping", {
  res <- suppressMessages(suppressWarnings(
    TESSERA_spNNGP(
      mock_sp_data, gene_name = "Gene1", em_iters = 20, 
      em_min_iters = 2, em_tol = 1.0, em_stopping = "abs_beta_norm",
      verbose = FALSE
    )
  ))
  
  expect_true(res$run_settings$em_iters_actual < 20)
})
