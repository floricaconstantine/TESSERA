## Functions to run after model fitting (inference).
# Dependencies in file: Matrix, dplyr, tibble.
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
summarizeTESSERAPerformance <- function(TESSERAData_obj, TESSERAOutput_obj) {
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
waldTestStastics <- function (TESSERAOutput_obj, contrast_mat) {
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

  V_hat <- inversePrecisionMatrixWald(TESSERAOutput_obj$beta_neghessian)
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
    # V_hat <- inversePrecisionMatrixWald(TESSERAOutput_obj$beta_neghessian[subset_idx, subset_idx])
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
inversePrecisionMatrixWald <- function(A) {
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


#' Fit a scaled non-central chi^2_1 distribution to a vector of statistics.
#'
#' We use the method of moments with the mean and variance to estimate parameters
#' of a scaled non-central chi squared distribution with one degree of freedom.
#' This fitting is NOT conditional.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param wald_stats Vector of Wald statistics: non-negative F-statistics.
#'  Square t-statistics to obtain F-statistics.
#' @param trimfrac NULL if ignored, c(LOWER, UPPER) if not.
#'  What fraction to trim from the left and right tails, e.g., c(0, 1e-4)
#'  clips nothing from the left tail and the top 1e-4 fraction from the right tail,
#'  i.e., everything above the 1 - 1e-4 quantile is trimmed.
#'
#' @returns Vector of (scaling, shift).
#'  (scaling) chi_1^2(shift).
#'
#' @note Consider using the trimfrac parameter if the data have outliers as
#'  this function is sensitive to the presence of large values/outliers.
#'
#' @importFrom stats var
#' @importFrom stats quantile
#' @export
#'
#' @examples
#' fit_scaled_noncentral_chi2(5 * stats::rchisq(10000, 1, ncp=10))
fit_scaled_noncentral_chi2 <- function(wald_stats, trimfrac = NULL) {
  # Toss out nan and inf
  wald_stats <- wald_stats[!is.nan(wald_stats)]
  wald_stats <- wald_stats[!is.infinite(wald_stats)]

  # Trim the vector of statistics
  if (!is.null(trimfrac)) {
    q_l <- stats::quantile(wald_stats, trimfrac[1])
    q_u <- stats::quantile(wald_stats, 1 - trimfrac[2])
    wald_stats <- wald_stats[wald_stats >= q_l]
    wald_stats <- wald_stats[wald_stats <= q_u]
  }

  v_w <- stats::var(wald_stats, na.rm = TRUE)
  m_w <- mean(wald_stats, na.rm = TRUE)

  scaling <- (1.0 / 2.0) * (2 * m_w - sqrt(2.0) * sqrt(max(0, 2.0 * m_w^2 - v_w)))
  shift <- sqrt(2.0) * sqrt(max(0, 2.0 * m_w^2 - v_w)) / (2.0 * scaling)

  return(c(scaling, shift))
}


#' Fit a scaled non-central chi^2_1 distribution to a vector of statistics.
#' The fit is conditional, below a threshold.
#'
#' We numerically maximize the likelihood.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param wald_stats Vector of Wald statistics: non-negative F-statistics.
#'  Square t-statistics to obtain F-statistics.
#' @param wald_thresh Threshold below which we fit.
#'
#' @returns Vector of (scaling, shift).
#'  (scaling) chi_1^2(shift).
#'
#' @note Consider using the trimfrac parameter if the data have outliers as
#'  this function is sensitive to the presence of large values/outliers.
#'
#' @importFrom DEoptim DEoptim DEoptim.control
#' @importFrom stats dchisq
#' @export
#'
#' @examples
#' fit_scaled_noncentral_chi2_DEoptim(5 * stats::rchisq(10000, 1, ncp=10), 20)
fit_scaled_noncentral_chi2_DEoptim <- function(wald_stats, wald_thresh) {
  # 1. Clean and subset data as before
  wald_stats <- wald_stats[is.finite(wald_stats)]
  wald_stats <- wald_stats[wald_stats < wald_thresh]

  if (length(wald_stats) < 10) {
    warning("Insufficient data points below threshold.")
    return(c(scaling = NA, shift = NA))
  }

  # 2. Internal Truncated NLL (Objective Function)
  # DEoptim expects the first argument to be the parameter vector
  trunc_nll <- function(params, x, thresh) {
    scale <- params[1]
    ncp <- params[2]

    # Protective boundary checks (DEoptim handles bounds, but this is a safety net)
    if (scale <= 1e-10 || ncp < 0)
      return(1e10)

    x[x == 0] <- 1e-8

    # Log-density calculation
    dens <- stats::dchisq(x / scale,
                          df = 1,
                          ncp = ncp,
                          log = TRUE) - log(scale)

    # Truncation correction (Log of the CDF up to threshold)
    log_F_thresh <- stats::pchisq(thresh / scale,
                                  df = 1,
                                  ncp = ncp,
                                  log.p = TRUE)

    val <- -(sum(dens) - length(x) * log_F_thresh)

    return(if (is.finite(val))
      val
      else
        1e10)
  }

  # 3. Define Global Search Bounds
  # Using your data characteristics: Scale ~670, NCP ~10-50
  lower_bounds <- c(1e-10, 0)
  upper_bounds <- c(max(1000, max(wald_stats) * 5), max(100, max(wald_stats)))

  # 4. Global Optimization
  fit_result <- tryCatch({
    DEoptim::DEoptim(
      fn = trunc_nll,
      lower = lower_bounds,
      upper = upper_bounds,
      x = wald_stats,
      thresh = wald_thresh,
      # DEoptim.control settings:
      # NP: Population size (usually 10 * length(params))
      # itermax: Number of generations
      control = DEoptim::DEoptim.control(
        NP = 20,
        itermax = 200,
        trace = FALSE # Set to TRUE if you want to see the "evolution"
      )
    )
  }, error = function(e) {
    warning("DEoptim failed: ", e$message)
    return(NULL)
  })

  if (is.null(fit_result))
    return(c(scaling = NA, shift = NA))

  # 5. Return the 'bestmem' (best member of the population)
  res <- fit_result$optim$bestmem
  names(res) <- c("scaling", "shift")
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
#' To use this function, the Wald statistics as well as a threshold are required.
#' The threshold is a rough upper bound for the statistics coming from the null
#' distribution--it need not be exact, but it should be greater than most of the
#' null distributed statistics and smaller than most of the non-null distributed
#' statistics.
#' In the absence of oracle knowledge, we recommend the `selectWaldStatisticThreshold`
#' function in this package for choosing an optimal threshold value.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param wald_stats Vector of Wald statistics: non-negative F-statistics.
#'  Square t-statistics to obtain F-statistics.
#' @param threshold A rough guess or upper bound for the Wald statistics under
#'  the null hypothesis.
#' @param conditional TRUE: p-values are conditional on X | X < Threshold;
#'  FALSE: p-values are not conditional.
#'
#' @returns Vector of p-values (un-adjusted for multiple corrections).
#'
#' @importFrom stats pchisq
#' @export
#'
#' @examples
#' waldStatisticPValuesThreshold(c(
#'   5 * stats::rchisq(2000, 1, ncp = 10),
#'   200 + 5 * stats::rchisq(2000, 1, ncp = 10)
#'  ), 200)
waldStatisticPValuesThreshold <- function(wald_stats, threshold, conditional =
                                            TRUE) {
  # Fit scaled, non-central chi-square distribution with one degree of freedom
  if (!conditional) {
    chi2_params <- fit_scaled_noncentral_chi2(wald_stats[wald_stats < threshold])
  } else {
    chi2_params <- fit_scaled_noncentral_chi2_DEoptim(wald_stats, threshold)
  }

  # RETURN UNCONDITIONAL P-VALUES
  return(stats::pchisq(
    wald_stats / chi2_params[1],
    1,
    ncp = chi2_params[2],
    lower.tail = FALSE
  ))

  ## THIS CODE IS NEVER RUN (IGNORE)
  # if (!conditional) {
  #   return(stats::pchisq(
  #     wald_stats / chi2_params[1],
  #     1,
  #     ncp = chi2_params[2],
  #     lower.tail = FALSE
  #   ))
  # } else {
  #   # Event A: P{X > observed}
  #   # Event B: P{X < threshold}
  #   # A intersect B: P{observed < X < threshold}
  #   #   = P{X < threshold} - P{X < observed}
  #   # P{A | B} = P{A intersect B} / P{B}
  #
  #   # P{B} = P{X < threshold}: Need to scale by chi2_params[1], as in
  #   # Wald statistics to match scaling
  #   p_B <- stats::pchisq(threshold / chi2_params[1],
  #                        1,
  #                        ncp = chi2_params[2],
  #                        lower.tail = TRUE)
  #
  #   # P{A intersect B} = P{observed < X < threshold}
  #   # First compute P{X < observed}
  #   p_A_B <- stats::pchisq(wald_stats / chi2_params[1],
  #                          1,
  #                          ncp = chi2_params[2],
  #                          lower.tail = TRUE)
  #   # Now take difference
  #   p_A_B <- p_B - p_A_B
  #   # Just in case, never needed
  #   p_A_B <- pmax(0, p_A_B)
  #   # Observed > threshold has probability zero under conditional
  #   p_A_B[wald_stats > threshold] <- 0
  #
  #   return(pmin(p_A_B / p_B, 1.0))
  # }
}

#' Given Wald statistics, compute an optimal threshold.
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
#' To use this function, the Wald statistics as well as a threshold is estimated.
#' The threshold is a rough upper bound for the statistics coming from the null
#' distribution--it need not be exact, but it should be greater than most of the
#' null distributed statistics and smaller than most of the non-null distributed
#' statistics.
#' The threshold is estimated by noting that under the null hypothesis,
#' the p-values should be approximately Uniform(0, 1)
#' distributed; for each threshold value, we compute the mean-squared error
#' between the quantiles of the p-values and those of the Uniform(0, 1) distribution
#' and choose the threshold with the lowest mean-squared error.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param wald_stats Vector of Wald statistics: non-negative F-statistics.
#'  Square t-statistics to obtain F-statistics.
#' @param trimfrac NULL if ignored, c(LOWER, UPPER) if not.
#'  What fraction to trim from the left and right tails, e.g., c(0, 1e-4)
#'  clips nothing from the left tail and the top 1e-4 fraction from the right tail,
#'  i.e., everything above the 1 - 1e-4 quantile is trimmed.
#'  In general, for this application, we recommend setting the lower/left value
#'  to 0 and the right value to something small to drop outliers.
#' @param conditional TRUE: p-values are conditional on X | X < Threshold;
#'  FALSE: p-values are not conditional.
#'
#' @returns Threshold value.
#'
#' @note Consider using the trimfrac parameter if the data have outliers:
#'  This function is sensitive to the presence of large values/outliers.
#'
#' @importFrom stats pchisq
#' @importFrom stats ppoints
#' @importFrom stats optimize
#' @importFrom stats quantile
#'
#' @export
#'
#' @examples
#' selectWaldStatisticThreshold(c(5 * stats::rchisq(2000, 1, ncp = 10),
#'  200 + 5 * stats::rchisq(2000, 1, ncp = 10)))
#' selectWaldStatisticThreshold(5 * stats::rchisq(2000, 1, ncp = 10))
selectWaldStatisticThreshold <- function(wald_stats,
                                         trimfrac = NULL,
                                         conditional = TRUE) {
  # Toss out nan and inf
  wald_stats <- wald_stats[!is.nan(wald_stats)]
  wald_stats <- wald_stats[!is.infinite(wald_stats)]

  # Trim the vector of statistics
  if (!is.null(trimfrac)) {
    q_l <- stats::quantile(wald_stats, trimfrac[1])
    q_u <- stats::quantile(wald_stats, 1 - trimfrac[2])
    wald_stats <- wald_stats[wald_stats >= q_l]
    wald_stats <- wald_stats[wald_stats <= q_u]
  }

  wrapper_optimization <- function(threshold) {
    # Subset to Wald statistics under threshold
    wald_subset <- wald_stats[wald_stats < threshold]

    # Guard against too-few statistics to actually fit anything
    if (length(wald_subset) < 5)
      return(1e10)

    # Calculate p-values
    wald_subset_p <- waldStatisticPValuesThreshold(wald_subset, threshold, conditional)

    # Obtain quantiles of the p-values
    # p_quantile <- stats::quantile(wald_subset_p, seq(0.0, 1.0, 1.0 / length(wald_subset_p)))

    theoretical <- stats::ppoints(length(wald_subset_p))
    mse <- mean((sort(wald_subset_p) - theoretical)^2)

    # Log MSE
    eps <- 1e-10
    mse_log <- mean((-log10(sort(wald_subset_p) + eps) --log10(theoretical + eps))^2)

    # MSE between observed quantiles and the Uniform[0, 1] distribution quantiles
    return(mse_log)
  }

  optim_out <- stats::optimize(wrapper_optimization, c(min(wald_stats), max(wald_stats)))
  return(optim_out$minimum)
}
