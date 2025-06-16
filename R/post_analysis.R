## Functions to run after model fitting (inference).
# Dependencies in file: Matrix, dplyr, tibble, moranfast.


#' Summarizes the results of the poisECM algorithms for a single gene.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param poisECMData_obj Object containing data.
#'  Created by the prepData method.
#'  Input to the poisECM algorithms.
#' @param poisECMOutput_obj Output of the poisECM algorithms.
#'
#' @returns A dataframe with MSE, spatial parameters, and Moran's I values.
#'  Moran's I values are only returned if the package moranfast is installed.
#'  Key (General):
#'    Each row corresponds to one sample.
#'    n_cells is the number of measurements in the sample.
#'    gamma_hat/tau2_hat are spatial parameters in lattice models.
#'    kernel_type, nugget_hat, range_hat, and smoothness_hat are spatial
#'      parameters in th spNNGP models.
#'  Key for MSE/means:
#'    _sample indicates that mean/MSE is computed for a single sample.
#'    _total indicates that mean/MSE is computed across all samples.
#'    counts2 refers to squared counts.
#'  Key for Moran's I values:
#'    counts: the observed counts/inputs.
#'    predictions: the predicted counts, equal to theta * (library size).
#'    residuals: counts - predictions.
#'    phi: the estimated spatial random effects.
#'    Xbeta: Covariates X multiplied by estimated effects beta.
#'    eta: phi + X beta.
#'    theta: Posterior expectation of exp(eta).
#'    librarysize: Often, the total counts per each cell.
#'
#' @note Requires the moranfast library if we want to return Moran's I of residuals.
#'    Install with `devtools::install_github('mcooper/moranfast')`.
#'
#' @importFrom tibble add_column
#' @importFrom dplyr bind_rows
#'
#' @export
summarizePoisECMPerformance <- function(poisECMData_obj, poisECMOutput_obj) {
  if (!requireNamespace("moranfast", quietly = TRUE)) {
    warning(
      "moranfast is not installed: This package requires the moranfast package to compute Moran's I values. Install it with devtools::install_github('mcooper/moranfast')"
    )
  }

  # Loop over samples to aggregate metrics
  inner_MSE_df <- list()
  for (s_idx in 1:length(poisECMData_obj$coords_list)) {
    # How to index the various vectors
    start_idx <- poisECMOutput_obj$start_idx_list[s_idx]
    if (s_idx < length(poisECMData_obj$coords_list)) {
      end_idx <- poisECMOutput_obj$start_idx_list[1 + s_idx] - 1
    } else {
      end_idx <- length(poisECMOutput_obj$phi_hat)
    }
    start_idx <- as.numeric(start_idx)
    end_idx <- as.numeric(end_idx)

    # Extract the relevant vectors
    # Predictions, residuals, random effects, X beta + phi (eta), X beta, theta
    counts <- poisECMData_obj$counts_list[[s_idx]][poisECMOutput_obj$run_settings$gene_idx, ]
    preds <- poisECMOutput_obj$predictions[start_idx:end_idx]
    resids <- poisECMOutput_obj$residuals[start_idx:end_idx]
    phis <- poisECMOutput_obj$phi_hat[start_idx:end_idx]
    etas <- poisECMOutput_obj$eta_hat[start_idx:end_idx]
    Xbs <- etas - phis
    thetas <- poisECMOutput_obj$theta_hat[start_idx:end_idx]
    library_sizes <- poisECMData_obj$library_size_list[[s_idx]]

    fit_model <- poisECMOutput_obj$run_settings$model_type
    if ("spNNGP" != fit_model) {
      gamma_hat <- poisECMOutput_obj$gamma_hat[s_idx]
      tau2_hat <- poisECMOutput_obj$tau2_hat[s_idx]
      kernel_type <- NA
      nugget_hat <- NA
      sill_hat <- NA
      range_hat <- NA
      smoothness_hat <- NA
    } else {
      if ("spNNGP" == fit_model) {
        gamma_hat <- NA
        tau2_hat <- NA
        kernel_type <- poisECMOutput_obj$run_settings$cov_type
        nugget_hat <- poisECMOutput_obj$cov_param_hat[s_idx, 1]
        sill_hat <- poisECMOutput_obj$cov_param_hat[s_idx, 2]
        range_hat <- poisECMOutput_obj$cov_param_hat[s_idx, 3]
        smoothness_hat <- poisECMOutput_obj$cov_param_hat[s_idx, 4]
      }
    }
    # Initialize the Moran's I vector
    # Only compute it if we have moranfast AND coordinates
    MI_vec <- rep(NA, 8)
    if (requireNamespace("moranfast", quietly = TRUE)) {
      if (!(
        is.null(poisECMData_obj$coords_list[[s_idx]][, 1]) ||
        is.null(poisECMData_obj$coords_list[[s_idx]][, 2])
      )) {
        MI_vec <- c(
          moranfast::calc_moran(
            counts,
            poisECMData_obj$coords_list[[s_idx]][, 1],
            poisECMData_obj$coords_list[[s_idx]][, 2]
          )[1],
          moranfast::calc_moran(
            preds,
            poisECMData_obj$coords_list[[s_idx]][, 1],
            poisECMData_obj$coords_list[[s_idx]][, 2]
          )[1],
          moranfast::calc_moran(
            resids,
            poisECMData_obj$coords_list[[s_idx]][, 1],
            poisECMData_obj$coords_list[[s_idx]][, 2]
          )[1],
          moranfast::calc_moran(
            phis,
            poisECMData_obj$coords_list[[s_idx]][, 1],
            poisECMData_obj$coords_list[[s_idx]][, 2]
          )[1],
          moranfast::calc_moran(
            etas,
            poisECMData_obj$coords_list[[s_idx]][, 1],
            poisECMData_obj$coords_list[[s_idx]][, 2]
          )[1],
          moranfast::calc_moran(
            Xbs,
            poisECMData_obj$coords_list[[s_idx]][, 1],
            poisECMData_obj$coords_list[[s_idx]][, 2]
          )[1],
          moranfast::calc_moran(
            thetas,
            poisECMData_obj$coords_list[[s_idx]][, 1],
            poisECMData_obj$coords_list[[s_idx]][, 2]
          )[1],
          moranfast::calc_moran(
            library_sizes,
            poisECMData_obj$coords_list[[s_idx]][, 1],
            poisECMData_obj$coords_list[[s_idx]][, 2]
          )[1]
        )
      }
    }
    inner_MSE_df[[s_idx]] <- data.frame(
      # Basic parameters
      gene = poisECMOutput_obj$run_settings$gene_name,
      fit_model = fit_model,
      sample = names(poisECMData_obj$coords_list)[s_idx],
      n_cells = (end_idx - start_idx + 1),
      # Spatial parameters
      gamma_hat = gamma_hat,
      tau2_hat = tau2_hat,
      kernel_type = kernel_type,
      nugget_hat = nugget_hat,
      range_hat = range_hat,
      smoothness_hat = smoothness_hat,
      # MSE of counts
      MSE_counts_sample = mean(resids^2),
      Mean_counts2_sample = mean(counts^2),
      # Moran's I of various quantities
      Moran_counts = MI_vec[1],
      Moran_predictions = MI_vec[2],
      Moran_residuals = MI_vec[3],
      Moran_phi = MI_vec[4],
      Moran_eta = MI_vec[5],
      Moran_Xbeta = MI_vec[6],
      Moran_theta = MI_vec[7],
      Moran_librarysize = MI_vec[8]
    )
  }
  inner_MSE_df <- dplyr::bind_rows(inner_MSE_df)
  # Add a column for the MSE across all samples
  inner_MSE_df <- tibble::add_column(
    inner_MSE_df,
    MSE_counts_total = sum(inner_MSE_df$MSE_counts_sample * inner_MSE_df$n_cells) / sum(inner_MSE_df$n_cells),
    .after = "MSE_counts_sample"
  )
  # Add a column for the mean counts across all samples
  inner_MSE_df <- tibble::add_column(
    inner_MSE_df,
    Mean_counts2_total = sum(inner_MSE_df$Mean_counts2_sample * inner_MSE_df$n_cells) / sum(inner_MSE_df$n_cells),
    .after = "Mean_counts2_sample"
  )

  return(inner_MSE_df)
}


#' Given the output of the poisECM algorithms and a contrast matrix, compute
#'  Wald T-statistics.
#'
#' @param poisECMOutput_obj Output of the poisECM algorithms.
#' @param contrast_mat Matrix with contrasts of the estimated coeficients.
#'  Rows are contrasts, columns correspond to columns in beta_hat (covariates).
#'
#' @returns A dataframe of Wald t-statistics.
#'
#' @importFrom dplyr bind_rows
#' @importFrom Matrix solve
#'
#' @export
waldTestStastics <- function (poisECMOutput_obj, contrast_mat) {
  # Extract coefficients and negative hessian and check dimensions
  beta_hat <- poisECMOutput_obj$beta_hat
  beta_names <- names(beta_hat)
  beta_hat_prec <- poisECMOutput_obj$beta_neghessian
  stopifnot(ncol(contrast_mat) == length(beta_hat))

  # Some overarching metadata
  fit_model <- poisECMOutput_obj$run_settings$model_type
  if ("spNNGP" != fit_model) {
    kernel_type <- NA
  } else {
    if ("spNNGP" == fit_model) {
      kernel_type <- poisECMOutput_obj$run_settings$cov_type
      if (is.null(kernel_type)) {
        kernel_type <- NA
      }
    } else {
      if (is.null(fit_model)) {
        fit_model <- NA
        kernel_type <- NA
      }
    }
  }
  gene <- poisECMOutput_obj$run_settings$gene_name
  if (is.null(gene)) {
    gene <- NA
  }

  wald_contrast_df <- list()
  for (c_idx in 1:nrow(contrast_mat)) {
    # Find indices involved in contrast and subset
    subset_idx <- which(contrast_mat[c_idx, ] != 0)
    R <- contrast_mat[c_idx, subset_idx]
    # Compute contrast value
    Rbeta <- sum(R * beta_hat[subset_idx])
    # Contrast SE
    RVR_inv <- sqrt(as.numeric(R %*% Matrix::solve(poisECMOutput_obj$beta_neghessian[subset_idx, subset_idx]) %*% R))

    wald_contrast_df[[1 + length(wald_contrast_df)]] <-
      data.frame(
        gene = gene,
        fit_model = fit_model,
        kernel_type = kernel_type,
        contrast_string = paste(R, beta_names[subset_idx], sep = "*", collapse = "+"),
        contrast_indices = paste(subset_idx, collapse = "_"),
        contrast_val = Rbeta,
        contrast_se = RVR_inv,
        wald_stat_t = Rbeta / RVR_inv
      )
  }

  wald_contrast_df <- dplyr::bind_rows(wald_contrast_df)
  rownames(wald_contrast_df) <- rownames(contrast_mat)
  return(wald_contrast_df)
}

