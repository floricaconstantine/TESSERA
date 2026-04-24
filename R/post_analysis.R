## Functions to run after model fitting (inference).
# Dependencies in file: Rcpp, Matrix, dplyr, tibble, optimx, pracma.
# Rcpp dependencies: calc_moran.cpp.


#' Summarizes the results of the TESSERA algorithms for a single gene.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param TESSERAData_obj Object containing data.
#'  Created by the prepData method.
#'  Input to the TESSERA algorithms.
#' @param TESSERAOutput_obj Output of the TESSERA algorithms.
#'
#' @returns A dataframe with MSE, spatial parameters, and Moran's I values.
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
#' @importFrom tibble add_column
#' @importFrom dplyr bind_rows
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @useDynLib TESSERA
#'
#' @export
summarize_TESSERA <- function(TESSERAData_obj, TESSERAOutput_obj) {
  # Loop over samples to aggregate metrics
  inner_MSE_df <- list()
  for (s_idx in 1:length(TESSERAData_obj$coords_list)) {
    # How to index the various vectors
    start_idx <- TESSERAOutput_obj$start_idx_list[s_idx]
    if (s_idx < length(TESSERAData_obj$coords_list)) {
      end_idx <- TESSERAOutput_obj$start_idx_list[1 + s_idx] - 1
    } else {
      end_idx <- length(TESSERAOutput_obj$phi_hat)
    }
    start_idx <- as.numeric(start_idx)
    end_idx <- as.numeric(end_idx)
    
    # Extract the relevant vectors
    # Predictions, residuals, random effects, X beta + phi (eta), X beta, theta
    counts <- TESSERAData_obj$counts_list[[s_idx]][TESSERAOutput_obj$run_settings$gene_idx, ]
    preds <- TESSERAOutput_obj$predictions[start_idx:end_idx]
    resids <- TESSERAOutput_obj$residuals[start_idx:end_idx]
    phis <- TESSERAOutput_obj$phi_hat[start_idx:end_idx]
    etas <- TESSERAOutput_obj$eta_hat[start_idx:end_idx]
    Xbs <- etas - phis
    thetas <- TESSERAOutput_obj$theta_hat[start_idx:end_idx]
    library_sizes <- TESSERAData_obj$library_size_list[[s_idx]]
    
    fit_model <- TESSERAOutput_obj$run_settings$model_type
    if ("spNNGP" != fit_model) {
      gamma_hat <- TESSERAOutput_obj$gamma_hat[s_idx]
      tau2_hat <- TESSERAOutput_obj$tau2_hat[s_idx]
      kernel_type <- NA
      nugget_hat <- NA
      sill_hat <- NA
      range_hat <- NA
      smoothness_hat <- NA
    } else {
      if ("spNNGP" == fit_model) {
        gamma_hat <- NA
        tau2_hat <- NA
        
        kernel_type <- TESSERAOutput_obj$run_settings$cov_type
        if (1 == length(TESSERAData_obj$coords_list)) {
          nugget_hat <- TESSERAOutput_obj$cov_param_hat[1]
          sill_hat <- TESSERAOutput_obj$cov_param_hat[2]
          range_hat <- TESSERAOutput_obj$cov_param_hat[3]
          smoothness_hat <- TESSERAOutput_obj$cov_param_hat[4]
        } else {
          nugget_hat <- TESSERAOutput_obj$cov_param_hat[s_idx, 1]
          sill_hat <- TESSERAOutput_obj$cov_param_hat[s_idx, 2]
          range_hat <- TESSERAOutput_obj$cov_param_hat[s_idx, 3]
          smoothness_hat <- TESSERAOutput_obj$cov_param_hat[s_idx, 4]
        }
      }
    }
    # Initialize the Moran's I vector
    # Only compute it if we have coordinates
    MI_vec <- rep(NA, 8)
    if (!(
      is.null(TESSERAData_obj$coords_list[[s_idx]][, 1]) ||
      is.null(TESSERAData_obj$coords_list[[s_idx]][, 2])
    )) {
      MI_vec <- c(
        calc_moran(
          counts,
          TESSERAData_obj$coords_list[[s_idx]][, 1],
          TESSERAData_obj$coords_list[[s_idx]][, 2]
        )[1],
        calc_moran(
          preds,
          TESSERAData_obj$coords_list[[s_idx]][, 1],
          TESSERAData_obj$coords_list[[s_idx]][, 2]
        )[1],
        calc_moran(
          resids,
          TESSERAData_obj$coords_list[[s_idx]][, 1],
          TESSERAData_obj$coords_list[[s_idx]][, 2]
        )[1],
        calc_moran(
          phis,
          TESSERAData_obj$coords_list[[s_idx]][, 1],
          TESSERAData_obj$coords_list[[s_idx]][, 2]
        )[1],
        calc_moran(
          etas,
          TESSERAData_obj$coords_list[[s_idx]][, 1],
          TESSERAData_obj$coords_list[[s_idx]][, 2]
        )[1],
        calc_moran(
          Xbs,
          TESSERAData_obj$coords_list[[s_idx]][, 1],
          TESSERAData_obj$coords_list[[s_idx]][, 2]
        )[1],
        calc_moran(
          thetas,
          TESSERAData_obj$coords_list[[s_idx]][, 1],
          TESSERAData_obj$coords_list[[s_idx]][, 2]
        )[1],
        calc_moran(
          library_sizes,
          TESSERAData_obj$coords_list[[s_idx]][, 1],
          TESSERAData_obj$coords_list[[s_idx]][, 2]
        )[1]
      )
    }
    
    inner_MSE_df[[s_idx]] <- data.frame(
      # Basic parameters
      gene = TESSERAOutput_obj$run_settings$gene_name,
      fit_model = fit_model,
      sample = names(TESSERAData_obj$coords_list)[s_idx],
      n_cells = (end_idx - start_idx + 1),
      # Spatial parameters
      gamma_hat = gamma_hat,
      tau2_hat = tau2_hat,
      kernel_type = kernel_type,
      nugget_hat = nugget_hat,
      sill_hat = sill_hat,
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


#' Given the output of the TESSERA algorithms and a contrast matrix, compute
#'  Wald T-statistics.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param TESSERAOutput_obj Output of the TESSERA algorithms.
#' @param contrast_mat Matrix with contrasts of the estimated coeficients.
#'  Rows are contrasts, columns correspond to columns in beta_hat (covariates).
#'  beta_hat is the vector of estimated coefficients, stored in TESSERAOutput_obj.
#'
#' @returns A dataframe of Wald t-statistics. Square these to get F-statistics.
#'
#' @importFrom dplyr bind_rows
#'
#' @export
#'
#' @examples
#' # Locate the saved model results in inst/extdata
#' rds_path <- system.file("extdata", "example_TESSERA_out_Leroux.rds", package = "TESSERA")
#' # Load the results object
#' TESSERA_out_Leroux <- readRDS(rds_path)
#'
#' # Run the Wald test
#' calc_Wald_statistics(TESSERAOutput_obj = TESSERA_out_Leroux,
#'                  contrast_mat = matrix(c(1, 0, 0), nrow = 1))
calc_Wald_statistics <- function (TESSERAOutput_obj, contrast_mat) {
  # Extract coefficients and negative hessian and check dimensions
  beta_hat <- TESSERAOutput_obj$beta_hat
  beta_names <- names(beta_hat)
  beta_hat_prec <- TESSERAOutput_obj$beta_neghessian
  stopifnot(ncol(contrast_mat) == length(beta_hat))
  
  # Some overarching metadata
  fit_model <- TESSERAOutput_obj$run_settings$model_type
  if ("spNNGP" != fit_model) {
    kernel_type <- NA
  } else {
    if ("spNNGP" == fit_model) {
      kernel_type <- TESSERAOutput_obj$run_settings$cov_type
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
  gene <- TESSERAOutput_obj$run_settings$gene_name
  if (is.null(gene)) {
    gene <- NA
  }
  
  V_hat <- invert_precision_matrix(TESSERAOutput_obj$beta_neghessian)
  wald_contrast_df <- list()
  for (c_idx in 1:nrow(contrast_mat)) {
    # Find indices involved in contrast and subset
    subset_idx <- which(contrast_mat[c_idx, ] != 0)
    # R <- contrast_mat[c_idx, subset_idx]
    R <- contrast_mat[c_idx, ]
    # Compute contrast value
    # Rbeta <- sum(R * beta_hat[subset_idx])
    Rbeta <- sum(R * beta_hat)
    
    # Contrast SE
    # V_hat <- invert_precision_matrix(TESSERAOutput_obj$beta_neghessian[subset_idx, subset_idx])
    RVR_inv <- sqrt(as.numeric(R %*% V_hat %*% R))
    
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


#' Helper function for Wald stat covariances.
#'
#' Invert a precision matrix, handling failure cases.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param A Estimated precision matrix.
#'
#' @returns Inverse of A.
#'
#' @import Matrix
#' @importFrom Matrix solve
#' @importFrom pracma pinv
invert_precision_matrix <- function(A) {
  # Try basic inversion first
  err_flag <- tryCatch({
    Ainv <- Matrix::solve(A)
    return(Ainv)
  }, error = function(cond) {
    warning("INVERSION OF Inv(V_hat) FAILED", "\n")
    err_flag <- TRUE
  })
  # Then try pseudoinverse
  if (err_flag) {
    err_flag <- tryCatch({
      Ainv <- pracma::pinv(as.matrix(A))
      return(Ainv)
    }, error = function(cond) {
      warning("PSEUDOINVERSE OF Inv(V_hat) FAILED", "\n")
      err_flag <- TRUE
    })
  }
  if (err_flag) {
    warning("Setting covariance to zero.", "\n")
    Ainv <- 0 * A
    return(Ainv)
  }
}


#' Fit a scaled non-central chi^2_1 distribution using BOBYQA
#'
#' @description
#' Numerically maximizes the likelihood of a scaled non-central chi-square
#' distribution (df=1) truncated above a specific threshold.
#' The fit is conditional on the threshold.
#'
#' @param wald_stats Vector of Wald statistics (non-negative).
#' @param wald_thresh Threshold below which data is included in the fit.
#'
#' @returns Vector of c(scale, shift).
#'
#' @importFrom optimx optimx
#' @importFrom stats dchisq pchisq quantile median
#' @export
#'
#' @examples
#' fit_scaled_noncentral_chi2(5 * stats::rchisq(1000, 1, ncp=10), 20)
fit_scaled_noncentral_chi2 <- function(wald_stats, wald_thresh) {
  # Data Cleaning
  wald_stats <- wald_stats[is.finite(wald_stats)]
  wald_stats <- wald_stats[wald_stats < wald_thresh]
  
  if (length(wald_stats) < 10) {
    warning("Insufficient data points below threshold.")
    return(c(scaling = NA, shift = NA))
  }
  
  # Truncated Negative Log-Likelihood (Objective Function)
  # optimx expects the first argument to be the parameter vector
  trunc_nll <- function(params, x, thresh) {
    scale <- params[1]
    ncp   <- params[2]
    
    # Safety check for log calculations
    x[x <= 0] <- 1e-8
    
    # Log-density: f(x; scale, ncp) = (1/scale) * f_chi2(x/scale; ncp)
    # log(f) = log_chi2_dens - log(scale)
    dens <- stats::dchisq(x / scale,
                          df = 1,
                          ncp = ncp,
                          log = TRUE) - log(scale)
    
    # Truncation correction: log(P(X < thresh))
    log_F_thresh <- stats::pchisq(thresh / scale,
                                  df = 1,
                                  ncp = ncp,
                                  log.p = TRUE)
    
    # Negative log-likelihood
    val <- -(sum(dens) - length(x) * log_F_thresh)
    
    return(if (is.finite(val))
      val
      else
        1e10)
  }
  
  # Method of Moments (MoM) Initialization
  v_w <- stats::var(wald_stats)
  m_w <- mean(wald_stats)
  
  # Your MoM logic
  term <- sqrt(max(0, 2.0 * m_w^2 - v_w))
  mom_scaling <- (2 * m_w - sqrt(2.0) * term) / 2.0
  mom_shift   <- (sqrt(2.0) * term) / (2.0 * max(1e-5, mom_scaling))
  
  # Define Bounds and Starting Values
  # Scale: 90% quantile of chi_1^2 is ~2.7.
  # We estimate upper scale by comparing observed 90th percentile to theoretical.
  upper_scale <- stats::quantile(wald_stats, 0.9) * (3.0 / 2.7)
  upper_ncp   <- min(10, stats::median(wald_stats))
  
  lower_bounds <- c(0.01, 0)
  upper_bounds <- c(max(2, upper_scale), max(2, upper_ncp))
  
  # Finalize Starting Parameters
  # Ensure MoM estimates aren't outside the box constraints
  start_params <- c(
    scaling = pmin(pmax(mom_scaling, lower_bounds[1]), upper_bounds[1]),
    shift   = pmin(pmax(mom_shift, lower_bounds[2]), upper_bounds[2])
  )
  
  # Optimization via optimx (bobyqa)
  fit_result <- tryCatch({
    optimx::optimx(
      par     = start_params,
      fn      = trunc_nll,
      lower   = lower_bounds,
      upper   = upper_bounds,
      method  = "bobyqa",
      x       = wald_stats,
      thresh  = wald_thresh,
      control = list(dowarn = FALSE)
    )
  }, error = function(e) {
    warning("optimx (bobyqa) failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(fit_result))
    return(c(scaling = NA, shift = NA))
  
  # Extract results
  # optimx returns a data frame; we want the parameters from the first row
  res <- c(scaling = fit_result$scaling[1], shift   = fit_result$shift[1])
  
  return(res)
}


#' Given Wald statistics, compute p-values.
#'
#' In the TESSERA algorithm, and in generalized linear mixed models in general,
#' the exact distribution of p-values for Wald statistics under the null hypothesis
#' is unknown.
#' We know that the standard errors are, for finite samples, likely biased
#' downward (too small).
#' Moreover, for finite samples, due to the approximations
#' in the E-step of the fitting algorithm, there may be a small positive bias.
#' Hence, a scaled chi-squared distribution or scaled F-distribution with
#' a non-centrality parameter might be a good model for the statistics under
#' the null hypothesis, where we note that an unscaled chi-squared or F distribution
#' are the standard asymptotic distributions for Wald statistics.
#' As we test a single hypothesis, we use a single degree of freedom.
#' To use this function, the Wald statistics as well as a threshold are required.
#' The threshold is a rough upper bound for the statistics coming from the null
#' distribution--it need not be exact, but it should be greater than most of the
#' null distributed statistics and smaller than most of the non-null distributed
#' statistics.
#' In the absence of oracle knowledge, we recommend the `select_Wald_threshold`
#' function in this package for choosing an optimal threshold value.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param wald_stats Vector of Wald statistics: non-negative F-statistics.
#'  Square t-statistics to obtain F-statistics.
#' @param threshold A rough guess or upper bound for the Wald statistics under
#'  the null hypothesis.
#'
#' @returns Vector of p-values (un-adjusted for multiple corrections).
#'
#' @importFrom stats pchisq
#' @export
#'
#' @examples
#' calc_Wald_pvalue_from_threshold(c(
#'   5 * stats::rchisq(2000, 1, ncp = 10),
#'   200 + 5 * stats::rchisq(2000, 1, ncp = 10)
#'  ), 200)
calc_Wald_pvalue_from_threshold <- function(wald_stats, threshold) {
  chi2_params <- fit_scaled_noncentral_chi2(wald_stats[wald_stats < threshold], threshold)
  
  # RETURN UNCONDITIONAL P-VALUES
  return(stats::pchisq(
    wald_stats / chi2_params[1],
    1,
    ncp = chi2_params[2],
    lower.tail = FALSE
  ))
}


#' Given Wald statistics, compute p-values.
#'
#' In the TESSERA algorithm, and in generalized linear mixed models in general,
#' the exact distribution of p-values for Wald statistics under the null hypothesis
#' is unknown.
#' We know that the standard errors are, for finite samples, likely biased
#' downward (too small).
#' Moreover, for finite samples, due to the approximations
#' in the E-step of the fitting algorithm, there may be a small positive bias.
#' Hence, a scaled chi-squared distribution or scaled F-distribution with
#' a non-centrality parameter might be a good model for the statistics under
#' the null hypothesis, where we note that an unscaled chi-squared or F distribution
#' are the standard asymptotic distributions for Wald statistics.
#' As we test a single hypothesis, we use a single degree of freedom.
#' To use this function, the Wald statistics as well as a threshold are required.
#' The threshold is a rough upper bound for the statistics coming from the null
#' distribution--it need not be exact, but it should be greater than most of the
#' null distributed statistics and smaller than most of the non-null distributed
#' statistics.
#' In the absence of oracle knowledge, we recommend the `select_Wald_threshold`
#' function in this package for choosing an optimal threshold value and parameters.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param wald_stats Vector of Wald statistics: non-negative F-statistics.
#'  Square t-statistics to obtain F-statistics.
#' @param chi2_params a vector of c(scale, shift) parameters.
#'
#' @returns Vector of p-values (un-adjusted for multiple corrections).
#'
#' @importFrom stats pchisq
#' @export
#'
#' @examples
#' calc_scaled_noncentral_chi2_pvalues(c(
#'   5 * stats::rchisq(2000, 1, ncp = 10),
#'   200 + 5 * stats::rchisq(2000, 1, ncp = 10)
#'  ), c(5, 10))
calc_scaled_noncentral_chi2_pvalues <- function(wald_stats, chi2_params) {
  # RETURN UNCONDITIONAL P-VALUES
  return(stats::pchisq(
    wald_stats / chi2_params[1],
    1,
    ncp = chi2_params[2],
    lower.tail = FALSE
  ))
}


#' Optimal threshold selection for empirical null estimation
#'
#' In generalized linear mixed models (GLMMs) and the TESSERA framework,
#' theoretical null distributions for Wald statistics are often unreliable
#' due to finite-sample biases and E-step approximations. This function
#' estimates an empirical null distribution—modeled as a scaled,
#' non-central \eqn{\chi^2_1} distribution—by identifying an optimal threshold
#' below which statistics likely originate from the null.
#'
#' The estimation procedure sweeps over a range of candidate thresholds. For
#' each threshold, it estimates parameters such that the resulting p-values
#' best approximate a \eqn{Uniform(0, 1)} distribution. The optimal threshold is
#' selected by minimizing a chosen error metric (e.g., MSE) between the
#' observed p-value quantiles and theoretical uniform quantiles.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param wald_stats A numeric vector of Wald statistics (non-negative).
#'   Note: t-statistics should be squared to obtain F-statistics (or 1-df \eqn{\chi^2}
#'   equivalent) before input.
#' @param quantile_spacing Numeric value defining the search resolution for
#'   the threshold. Defaults to 0.01 (searching 0.01, 0.02, ..., 0.99 quantiles).
#' @param metric Character string specifying the selection metric.
#'   Options include "Raw_" or "Log_" prefixed to "MSE", "MAE", "MedAE",
#'   or "MaxAE". These correspond to Mean Squared Error, Mean Absolute Error,
#'   Median Absolute Error, and Maximum Absolute Error (L-infinity), respectively.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \strong{threshold}: The selected optimal Wald statistic threshold.
#'   \item \strong{chi2params}: A vector containing the estimated shift (non-centrality)
#'     and scale parameters for the empirical null.
#'   \item \strong{threshold_results}: A data frame summarizing metrics and
#'     fitted parameters across all evaluated threshold quantiles.
#' }
#'
#' @importFrom stats pchisq ppoints quantile
#' @importFrom dplyr bind_rows
#'
#' @export
#'
#' @examples
#' # Simulate null and non-null Wald statistics
#' null_stats <- 5 * stats::rchisq(2000, df = 1, ncp = 0.5)
#' alt_stats <- 200 + 5 * stats::rchisq(500, df = 1, ncp = 10)
#' wald_vec <- c(null_stats, alt_stats)
#'
#' # Select optimal threshold
#' select_Wald_threshold(wald_vec, quantile_spacing = 0.05)
select_Wald_threshold <- function (wald_stats,
                                   quantile_spacing = 0.01,
                                   metric = "Raw_MSE") {
  # Function to get errors and fits at a given threshold value
  wrapper_function <- function(wald_thresh, wald_data) {
    # Subset to statistics under threshold
    wald_nulls <- wald_data[wald_data < wald_thresh]
    
    # Fit distribution (conditionally) to statistics under threshold
    chi2_fit <- fit_scaled_noncentral_chi2(wald_nulls, wald_thresh)
    # Just in case
    if (any(is.na(chi2_fit)))
      return(NULL)
    
    # Compute p-values
    # We use the fitted scaling (chi2_fit[1]) and ncp (chi2_fit[2])
    wald_null_pvals <- stats::pchisq(
      wald_nulls / chi2_fit[1],
      df = 1,
      ncp = chi2_fit[2],
      lower.tail = FALSE
    )
    
    # Prepare Quantiles for comparison
    obs <- sort(wald_null_pvals)
    theo <- stats::ppoints(length(obs))
    
    # Logged versions
    eps <- 1e-10 # Small constant to prevent log(0)
    obs_log  <- -log10(obs + eps)
    theo_log <- -log10(theo + eps)
    
    # Helper function to calculate error metrics
    calc_errors <- function(o, t) {
      err <- abs(o - t)
      list(
        MSE   = mean(err^2),
        MAE   = mean(err),
        MedAE = median(err),
        MaxAE = max(err)
      )
    }
    # Calculate for both scales
    raw_errors <- calc_errors(obs, theo)
    log_errors <- calc_errors(obs_log, theo_log)
    
    # Combine into a single dataframe
    res <- data.frame(
      threshold = wald_thresh,
      scale     = as.numeric(chi2_fit[1]),
      shift     = as.numeric(chi2_fit[2]),
      n_obs     = length(wald_nulls)
    )
    # Append metrics with prefixes
    res[, paste0("Raw_", names(raw_errors))] <- as.list(raw_errors)
    res[, paste0("Log_", names(log_errors))] <- as.list(log_errors)
    
    return(res)
  }
  
  # Quantiles to use for threshold
  quantile_list <- seq(0 + quantile_spacing, 1.0 - quantile_spacing, quantile_spacing)
  
  # Exclude NA statistics
  wald_stats <- wald_stats[!is.na(wald_stats)]
  threshold_grid  <- stats::quantile(wald_stats, quantile_list, na.rm = TRUE)
  
  # Loop across Thresholds
  fit_results_list <- list()
  for (t_idx in 1:length(threshold_grid)) {
    fit_results_list[[t_idx]] <- wrapper_function(threshold_grid[t_idx], wald_stats)
  }
  # Combine the list into one dataframe
  fit_results_df <- dplyr::bind_rows(fit_results_list)
  # Add in quantile
  fit_results_df$quantile <- quantile_list
  
  # Find which threshold minimizes the metric
  min_idx <- which.min(fit_results_df[, metric])
  
  return(
    list(
      threshold = fit_results_df$threshold[min_idx],
      chi2_params = c(scale = fit_results_df$scale[min_idx], shift =
                        fit_results_df$shift[min_idx]),
      threshold_results = fit_results_df
    )
  )
}
