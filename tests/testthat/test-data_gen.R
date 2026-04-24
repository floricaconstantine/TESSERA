library(testthat)
library(Matrix)

# --- Shared Setup ---
set.seed(2026)
beta_true <- c(1, -0.5)
n_pts <- 50
nb_dist <- 0.2

test_that("generate_data_one_area produces correct structure and dimensions",
          {
            # Test Leroux model with Bernoulli covariates
            res <- generate_data_one_area(
              n_points = n_pts,
              nb_dist = nb_dist,
              model_type = "Leroux",
              beta_true = beta_true,
              gamma_true = 0.5,
              tau2_true = 1.0,
              X_type = "rand_bern"
            )
            
            # Check required output fields
            expect_setequal(
              names(res),
              c(
                "X",
                "W",
                "D",
                "eig_list",
                "library_size",
                "x_coords",
                "y_coords",
                "coords",
                "Q",
                "phi_true",
                "eta_true",
                "theta_true",
                "z"
              )
            )
            
            # Check dimensions
            expect_equal(length(res$z), n_pts)
            expect_equal(nrow(res$X), n_pts)
            expect_equal(dim(res$W), c(n_pts, n_pts))
            
            # Verify coordinate sorting (X then Y)
            expect_false(is.unsorted(res$coords[, 1]))
            
            # Verify adjacency symmetry
            expect_true(isSymmetric(res$W))
          })

test_that("generate_data_one_area_spNNGP handles spatial GP data", {
  cov_params <- c(0.1, 1.0, 0.5) # nugget, sill, range
  
  res <- generate_data_one_area_spNNGP(
    n_points = n_pts,
    nb_dist = nb_dist,
    cov_type = "Exp",
    cov_params = cov_params,
    nngp_k = 5,
    beta_true = beta_true
  )
  
  expect_s4_class(res$Q, "Matrix")
  expect_equal(length(res$phi_true), n_pts)
  
  # FIX: Nest the expectations to catch BOTH warnings emitted by the function
  expect_warning(
    expect_warning(
      res_zero <- generate_data_one_area_spNNGP(
        n_points = 10,
        nb_dist = 0.5,
        cov_type = "Exp",
        cov_params = c(0, 0, 1),
        nngp_k = 3,
        beta_true = 1
      ),
      "Q IS ZERO",
      fixed = TRUE
    ),
    "PHI IS ZERO",
    fixed = TRUE
  )
  
  expect_true(all(res_zero$phi_true == 0))
})

test_that("sample_Poisson_lattice is consistent with generated structures",
          {
            # Generate a base structure
            base <- generate_data_one_area(20, 0.3, "CAR", 1, 0.5, 1, "intercept")
            
            # Sample new counts using existing structure
            res <- sample_Poisson_lattice(
              model_type = "CAR",
              X = base$X,
              W = base$W,
              D = base$D,
              library_size = base$library_size,
              tau2_true = 1.2,
              gamma_true = 0.6,
              beta_true = 1.5
            )
            
            expect_named(res, c("z", "phi_true", "eta_true", "theta_true"))
            expect_equal(length(res$z), 20)
          })

test_that("prep_synth_data integrates counts into TESSERA objects", {
  # Create a minimal mock TESSERAData object
  # (Assuming the object has counts_list, X_list, W_list, etc.)
  mock_obj <- list(
    counts_list = list(S1 = matrix(rpois(10, 5), nrow = 1)),
    X_list = list(S1 = matrix(1, 10, 1)),
    W_list = list(S1 = Diagonal(10)),
    D_list = list(S1 = Diagonal(10)),
    library_size_list = list(S1 = rep(1, 10))
  )
  # Mock names for genes
  rownames(mock_obj$counts_list$S1) <- "GeneA"
  
  expect_error(
    prep_synth_data(mock_obj, "GeneA", "INVALID_MODEL", beta_true = 1),
    "Invalid data_gen_model"
  )
})
