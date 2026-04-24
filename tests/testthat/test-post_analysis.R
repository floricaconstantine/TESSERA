library(testthat)
library(Matrix)
library(dplyr)
library(tibble)

# --- Shared Setup ---
set.seed(2026)

mock_data <- structure(list(
  coords_list = list(S1 = matrix(runif(20), ncol = 2)),
  counts_list = list(S1 = matrix(rpois(10, 5), nrow = 1)),
  library_size_list = list(S1 = rep(1, 10))
),
class = "TESSERAData")

mock_output <- structure(
  list(
    start_idx_list = c(S1 = 1),
    predictions = rpois(10, 5),
    residuals = rnorm(10),
    phi_hat = rnorm(10),
    eta_hat = rnorm(10),
    theta_hat = exp(rnorm(10)),
    beta_hat = c(intercept = 1.2, treatment = -0.5),
    beta_neghessian = matrix(c(10, 2, 2, 5), 2, 2),
    gamma_hat = 0.5,
    tau2_hat = 0.8,
    run_settings = list(
      gene_idx = 1,
      gene_name = "GeneA",
      model_type = "CAR"
    )
  ),
  class = "TESSERAOutput"
)
names(mock_output$beta_hat) <- c("intercept", "treatment")

# --- Tests ---

test_that("summarize_TESSERA aggregates metrics correctly", {
  # suppressMessages silences internal Moran's I logs
  res <- suppressMessages(summarize_TESSERA(mock_data, mock_output))
  
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_true("Moran_phi" %in% colnames(res))
})

test_that("calc_Wald_statistics computes correct dimensions", {
  contrast <- matrix(c(0, 1), nrow = 1)
  res <- calc_Wald_statistics(mock_output, contrast)
  
  expect_s3_class(res, "data.frame")
  # Use as.numeric to strip the 'treatment' name coming from beta_hat
  expect_equal(
    as.numeric(res$wald_stat_t),
    as.numeric(mock_output$beta_hat[2] / res$contrast_se)
  )
})

test_that("invert_precision_matrix handles singular matrices", {
  A_good <- matrix(c(2, 0, 0, 2), 2, 2)
  expect_silent(invert_precision_matrix(A_good))
  
  # When A is all zeros, solve() fails (Warning 1), then pinv() succeeds.
  # We catch the first warning that actually fires.
  A_bad <- matrix(0, 2, 2)
  expect_warning(
    invert_precision_matrix(A_bad),
    "INVERSION.*FAILED" # Regex match to handle potential newlines
  )
})

test_that("fit_scaled_noncentral_chi2 handles insufficient data", {
  # Providing 5 points (needs 10) triggers the fallback NA return[source: 6]
  expect_warning(res <- fit_scaled_noncentral_chi2(rnorm(5), 10),
                 "Insufficient data")
  expect_true(all(is.na(res)))
})

test_that("select_Wald_threshold identifies an optimal threshold", {
  # Large dataset to ensure the optimizer converges and doesn't spam NaNs
  null_stats <- 2 * stats::rchisq(200, df = 1, ncp = 2)
  signal_stats <- stats::rchisq(50, df = 1, ncp = 50)
  wald_all <- c(null_stats, signal_stats)
  
  # Silence numerical warnings from the optimizer during threshold search
  res <- suppressWarnings(select_Wald_threshold(wald_all, quantile_spacing = 0.2))
  
  expect_named(res, c("threshold", "chi2_params", "threshold_results"))
  expect_gt(res$threshold, 0)
})

test_that("calc_Wald_pvalue_from_threshold computes p-values properly",
          {
            # Stats vector must have >10 points below threshold to satisfy fit_scaled_noncentral_chi2[source: 6]
            stats <- c(stats::rchisq(15, 1, ncp = 2), 100, 200)
            pvals <- calc_Wald_pvalue_from_threshold(stats, threshold = 50)
            
            expect_length(pvals, 17)
            expect_true(all(pvals >= 0 & pvals <= 1, na.rm = TRUE))
          })
