## Wrapper functions to fit model.
# Dependencies in file: Matrix.
# Dependencies: Functions from utils.R, models.R, E_step.R, M_step.R.
# Dependences from functions in other files not listed: pracma, sp, gstat.
# Rcpp dependencies: calc_moran.cpp.


#' Fit multi-area Poisson spatial generalized linear model.
#'  Allows for sparse nearest neighbor Gaussian process random effects.
#'  Exponential, Gaussian, Matern, and Spherical covariance structures
#'  are supported.
#'  Fits a common set of fixed effects across areas.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param TESSERAData_obj Object containing data.
#'  Created by the prepData method.
#' @param gene_name Which gene/measurement to fit.
#'  I.e., Which row in the count data matrix to fit.
#' @param cov_type Which model to fit for the Gaussian process covariance.
#'    "Exp", "Mat", "Gau", and "Sph" are the valid options, for
#'    Exponential, Matern, Gaussian, and Spherical, respectively.
#' @param nngp_k How many neighbors to use for the spNNGP kernel.
#' @param em_iters Number of E(C)M iterations to run.
#' @param opt_iters Number of inner M-step iterations to run.
#'    The M-step is technically a CM step, i.e., a conditional maximization.
#'    That is, we perform block coordinate ascent wherein we maximize the
#'    expected likelihood in tau^2 and then in gamma for each area, and then
#'    in beta, where each maximization holds all other parameters constant.
#'    This approach leads to a slower convergence in terms of the number of
#'    iterations (em_iters) but faster computation in each iteration.
#' @param em_min_iters Minimum number of E(C)M iterations to run.
#' @param em_tol Tolerance for early stopping the E(C)M algorithm.
#' @param em_stopping How to decide to stop early.
#'    NULL: Don't stop early.
#'    "abs_loglike": Difference in absolute data log likelihood, |old - new|.
#'    "rel_loglike": Relative difference in data log likelihood, |old - new|/|old|.
#'      Data log likelihood taken across alll areas.
#'    "abs_beta_norm": Absolute L2 norm of difference in beta, ||old - new||_2.
#'    "rel_beta_norm": Relative L2 norm of difference in beta,
#'      ||old - new||_2 / ||old||_2.
#' @param beta_init How to initialize beta.
#'    "glm": Stack up all the observed counts z and covariate matrices X and
#'      fit a Poisson GLM z ~ 0 + X using base R's glm function.
#'    "random": Initialize with standard normal random variables.
#'    A vector or matrix: Initialize with provided value.
#' @param cov_init How to initialize covariance parameters.
#'    "BRISC": Call BRISC on log(0.5 + z) - X beta_init.
#'    "variogram": log(0.5 + z) - X beta_init.
#'    A vector of matrix: Initialize with provided value.
#' @param cov_fit_method "BRISC" or "variogram": How to perform the M-step
#'    estimation/update of covariance parameters.
#' @param verbose Boolean of whether to print out updates.
#' @param dense_matrices Boolean of whether to treat Q as a dense matrix when
#'    computing certain quantities in the E-step. This will increase memory usage
#'    (e.g., ~2 -> 8 GB for 3k cells per sample and 20 neighbors), but will dramatically
#'    speed up computation (around a 60-80% decrease).
#'
#' @return A list comprised of the following:
#' @return beta_hat: Estimated coefficients.
#' @return cov_param_hat: Estimated spatial covariance parameters.
#' @returns predictions: Predicted values.
#'    Vector of all values (corresponds to z_list stacked together).
#' @returns eta_hat: Estimated random effects.
#'    Vector of all values (corresponds to z_list stacked together).
#' @returns theta_hat: Estimated Poisson parameters.
#'    Vector of all values (corresponds to z_list stacked together).
#' @returns phi_hat: Estimated spatial random effects.
#'    Vector of all values (corresponds to z_list stacked together).
#' @returns cov_param_tracker: Estimated correlation parameters.
#'    areas x EM iterations x (Nugget, Sill, Range, Smoothness)---history across iterations.
#' @returns beta_tracker: Estimated fixed effects.
#'    coefficients x EM iterations matrix---history across iterations.
#' @returns R2_tracker: Squared correlation between predictions and observations.
#'    areas x EM Iterations matrix---history across iterations.
#' @returns MSE_tracker: Mean-Square Error between predictions and observations.
#'    areas x EM Iterations matrix---history across iterations.
#' @returns data_log_like_tracker: Log likelihood of observations.
#'    areas x EM Iterations matrix---history across iterations.
#' @returns expected_log_like_tracker: Expected log likelihood at current parameters.
#'    areas x EM Iterations matrix---history across iterations.
#' @returns resid_moran_nb: Moran's I for each area computed using neighbor
#'   adjacency as weights.
#'   areas x EM Iterations matrix---history across iterations.
#' @returns resid_moran: Moran's I for each area computed using coordinates.
#'   areas x EM Iterations x (value, expectation, sd) array.
#'   NULL unless coordinates are supplied.
#' @returns start_idx_list: Indices where each area's values start in predictions, etc.
#'    E.g., if area 1 has 100 points and area 2 has 50, we would have
#'    c(1, 101, 151, ...).
#' @returns time: Total time taken by function.
#' @returns beta_neghessian: Negative Hessian of final value of beta.
#'    Matrix.
#' @returns run_settings: A list with the parameter settings used to run the algorithm.
#' @returns performanceSummary: A dataframe with summary statistics for each sample.
#'  See [summarizeTESSERAPerformance()] for more details.
#'
#' @references Meng, Xiao-Li, and Donald B. Rubin.
#'                "Maximum likelihood estimation via the ECM algorithm: A general framework."
#'                Biometrika 80.2 (1993): 267-278.
#'
#' @note Requires the Matrix library.
#'
#' @import Matrix
#' @importFrom stats coef
#' @importFrom stats cor
#' @importFrom stats dpois
#' @importFrom stats poisson
#' @importFrom stats rnorm
#' @importFrom stats var
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @useDynLib TESSERA
#'
#' @export
#' 
#' @example inst/examples/gen_data_and_run.R
TESSERA_spNNGP <- function(TESSERAData_obj,
                           gene_name,
                           cov_type = "Exp",
                           nngp_k = 20,
                           em_iters = 200,
                           opt_iters = 5,
                           em_min_iters = 15,
                           em_tol = 1e-3,
                           em_stopping = NULL,
                           beta_init = "glm",
                           cov_init = "BRISC",
                           cov_fit_method = "BRISC",
                           verbose = FALSE,
                           dense_matrices = TRUE)
{
  # Start clock
  t0_EM <- Sys.time()
  
  # Check inputs
  checkInputsTESSERAspNNGP(TESSERAData_obj)
  
  # Extract gene of interest and associated counts
  gene_idx <- which(rownames(TESSERAData_obj$counts_list[[1]]) == gene_name)
  z_list <- lapply(TESSERAData_obj$counts_list, function (x) {
    x[gene_idx, ]
  })
  
  # Number of areas
  n_areas <- length(z_list)
  # Total number of points
  n_total_points <- 0
  # Starting points in larger matrices
  start_idx_list <- rep(NA, length(z_list))
  for (idx in 1:length(TESSERAData_obj$X_list)) {
    n_total_points <- n_total_points + nrow(TESSERAData_obj$X_list[[idx]])
    if (1 == idx) {
      start_idx_list[idx] <- 1
    } else {
      start_idx_list[idx] <- start_idx_list[idx - 1] + length(z_list[[idx - 1]])
    }
  }
  # Dimension of coefficient vector
  beta_dim <- ncol(TESSERAData_obj$X_list[[1]])
  
  # Set up nearest neighbor tracker and associated distances
  sp_dist_list <- list()
  for (idx in 1:length(TESSERAData_obj$coords_list)) {
    sp_dist_list[[idx]] <- sparseDist(TESSERAData_obj$coords_list[[idx]], nngp_k)
  }
  
  # Store parameter estimates and track
  cov_param_tracker <- array(data = NA, dim = c(n_areas, 1 + em_iters, 4))
  beta_tracker <- matrix(data = NA,
                         nrow = beta_dim,
                         ncol = em_iters + 1)
  fit_tracker <- matrix(data = NA,
                        nrow = n_total_points,
                        ncol = em_iters)
  eta_tracker <- matrix(data = NA,
                        nrow = n_total_points,
                        ncol = em_iters)
  theta_tracker <- matrix(data = NA,
                          nrow = n_total_points,
                          ncol = em_iters)
  
  # Track performance
  R2_tracker <- matrix(data = NA,
                       nrow = n_areas,
                       ncol = em_iters)
  MSE_tracker <- matrix(data = NA,
                        nrow = n_areas,
                        ncol = em_iters)
  data_log_like_tracker <- matrix(data = NA,
                                  nrow = n_areas,
                                  ncol = em_iters)
  expected_log_like_tracker <- matrix(data = NA,
                                      nrow = n_areas,
                                      ncol = em_iters)
  resid_moran <- array(data = NA, dim = c(n_areas, em_iters, 3))
  
  # Initialize parameters: beta
  if (is.character(beta_init) && ("random" == beta_init)) {
    beta_tracker[, 1] <- stats::rnorm(beta_dim)
    message("Random initialization for beta.", "\n")
  } else if (is.character(beta_init) && ("glm" == beta_init)) {
    # Reasonable initialization: a basic GLM
    
    # Stack into a single vector/matrix
    z_vec <- Reduce(c, z_list)
    lib_vec <- Reduce(c, TESSERAData_obj$library_size_list)
    X_mat <- Reduce(rbind, TESSERAData_obj$X_list)
    # Fit GLM
    beta_tmp <- as.vector(stats::coef(
      # stats::glm(z_vec / lib_vec ~ 0 + X_mat, family = stats::poisson())
      stats::glm(
        z_vec ~ 0 + X_mat,
        family = stats::poisson(),
        offset = log(lib_vec)
      )
    ))
    beta_tmp[is.nan(beta_tmp)] <- 0
    beta_tmp[is.infinite(beta_tmp)] <- 0
    beta_tracker[, 1] <- beta_tmp
    
    # Memory
    rm(z_vec)
    rm(X_mat)
    
    message("GLM initialization for beta.", "\n")
  } else if (is.numeric(beta_init) &&
             (is.vector(beta_init) || is.matrix(beta_init))) {
    # Pass in a value
    beta_tracker[, 1] <- as.vector(beta_init)
    
    message("Pre-defined initialization for beta.", "\n")
  } else {
    stop("Invalid initialization for beta.")
  }
  # Handle NA
  beta_tracker[is.nan(beta_tracker[, 1]), 1] <- 0.0
  beta_tracker[is.na(beta_tracker[, 1]), 1] <- 0.0
  message("Initial beta ", paste(beta_tracker[, 1], collapse = " "), "\n")
  
  # Initialize parameters: covariance structure
  if (is.character(cov_init) && ("BRISC" == cov_init)) {
    for (area_idx in 1:n_areas) {
      param_est <- M_step_BRISC(
        log((0.5 + z_list[[area_idx]]) / TESSERAData_obj$library_size_list[[area_idx]]),
        beta_tracker[, 1],
        TESSERAData_obj$X_list[[area_idx]],
        TESSERAData_obj$coords_list[[area_idx]],
        cov_type,
        nngp_k
      )
      cov_param_tracker[area_idx, 1, 1:length(param_est)] <- param_est
    }
    message("BRISC initialization for covariances.", "\n")
  } else if (is.character(cov_init) && ("variogram" == cov_init)) {
    for (area_idx in 1:n_areas) {
      param_est <- M_step_variogram(
        log((0.5 + z_list[[area_idx]]) / TESSERAData_obj$library_size_list[[area_idx]]),
        beta_tracker[, 1],
        TESSERAData_obj$X_list[[area_idx]],
        TESSERAData_obj$coords_list[[area_idx]],
        cov_type
      )
      cov_param_tracker[area_idx, 1, 1:length(param_est)] <- param_est
    }
    message("Variogram initialization for covariances.", "\n")
  } else if (is.vector(cov_init)) {
    for (area_idx in 1:n_areas) {
      cov_param_tracker[area_idx, 1, 1:length(cov_init)] <- cov_init
    }
    message("Predefined initialization for covariances: Common.", "\n")
  } else if (is.matrix(cov_init)) {
    cov_param_tracker[, 1, 1:ncol(cov_init)] <- cov_init
    message("Predefined initialization for covariances: Individual.",
            "\n")
  } else {
    stop("Invalid initialization for parameters.")
  }
  message("Initial covariance parameters ",
          paste(cov_param_tracker[, 1, ], collapse = " "),
          "\n")
  
  
  # Also initialize dependency Q as a function of the variogram parameters
  Q_hat_list <- list()
  Dinv_list <- list()
  A_hat_list <- list()
  for (area_idx in 1:n_areas) {
    param_est <- nngp_prec_mat(sp_dist_list[[area_idx]],
                               TESSERAData_obj$coords_list[[area_idx]],
                               cov_type,
                               cov_param_tracker[area_idx, 1, ])
    Q_hat_list[[area_idx]] <- param_est$Q
    Dinv_list[[area_idx]] <- param_est$Dinv
    A_hat_list[[area_idx]] <- param_est$A
  }
  
  # Loop over EM iterations
  for (em_idx in 1:em_iters) {
    if (verbose || (0 == (em_idx %% 100))) {
      message(
        paste(
          "Start of EM Iteration ",
          em_idx,
          " of ",
          em_iters,
          "; ",
          Sys.time() - t0_EM,
          " Elapsed"
        ),
        "\n"
      )
    }
    
    # Run E-Step: Get covariance and mean of eta
    Vhat_list <- list()
    eta_hat_list <- list()
    for (area_idx in 1:n_areas) {
      if (dense_matrices) {
        Vhat_list[[area_idx]] <- E_step_Vhat(as.matrix(Q_hat_list[[area_idx]]), 1.0, z_list[[area_idx]])
      } else {
        Vhat_list[[area_idx]] <- E_step_Vhat(Q_hat_list[[area_idx]], 1.0, z_list[[area_idx]])
      }
      
      eta_hat_list[[area_idx]] <- E_step_etahat(
        Vhat_list[[area_idx]],
        Q_hat_list[[area_idx]],
        1.0,
        beta_tracker[, em_idx],
        TESSERAData_obj$X_list[[area_idx]],
        z_list[[area_idx]],
        TESSERAData_obj$library_size_list[[area_idx]]
      )
      # Get a few more things out of the E-step: theta and predictions
      theta_hat <- as.numeric(E_step_thetahat(Vhat_list[[area_idx]], eta_hat_list[[area_idx]]))
      z_hat <- as.numeric(E_step_predict(theta_hat, TESSERAData_obj$library_size_list[[area_idx]]))
      
      # Store stuff
      fit_tracker[start_idx_list[area_idx]:(start_idx_list[area_idx] + length(z_hat) - 1), em_idx] <- z_hat
      eta_tracker[start_idx_list[area_idx]:(start_idx_list[area_idx] + length(eta_hat_list[[area_idx]]) - 1), em_idx] <- as.numeric(eta_hat_list[[area_idx]])
      theta_tracker[start_idx_list[area_idx]:(start_idx_list[area_idx] + length(theta_hat) - 1), em_idx] <- theta_hat
      
      # Performance
      R2_tracker[area_idx, em_idx] <- stats::cor(z_hat, z_list[[area_idx]])^2
      MSE_tracker[area_idx, em_idx] <- mean(abs(z_hat - z_list[[area_idx]])^2)
      # Observed, expected, sd
      resid_moran[area_idx, em_idx, ] <- calc_moran(
        z_hat - z_list[[area_idx]],
        TESSERAData_obj$coords_list[[area_idx]][, 1],
        TESSERAData_obj$coords_list[[area_idx]][, 2]
      )
      
      
      # Data log likelihood
      # data_log_like_tracker[area_idx, em_idx] <- poisson_loglike(z_list[[area_idx]], theta_hat * library_size_list[[area_idx]])
      data_log_like_tracker[area_idx, em_idx] <- sum(
        stats::dpois(
          round(z_list[[area_idx]]),
          theta_hat * TESSERAData_obj$library_size_list[[area_idx]],
          log = TRUE
        ),
        na.rm = TRUE
      )
      
      # Expected log likelihood
      expected_log_like_tracker[area_idx, em_idx] <- expected_loglike(
        Vhat_list[[area_idx]],
        eta_hat_list[[area_idx]],
        Q_hat_list[[area_idx]],
        0.0,
        # gamma
        1.0,
        # tau^2
        beta_tracker[, em_idx],
        TESSERAData_obj$X_list[[area_idx]],
        A_hat_list[[area_idx]],
        # W
        A_hat_list[[area_idx]],
        # D
        Dinv_list[[area_idx]],
        # Eigenvalues
        "spNNGP"
      )
    }
    
    # Run (C)M-Step: Optimize
    for (opt_idx in 1:opt_iters) {
      # Initialize with current values
      if (1 == opt_idx) {
        cov_param_tracker[area_idx, 1 + em_idx, ] <- cov_param_tracker[area_idx, em_idx, ]
        beta_tracker[, 1 + em_idx] <- beta_tracker[, em_idx]
      }
      
      for (area_idx in 1:n_areas) {
        # Optimize in covariance parameters
        if ("variogram" == cov_fit_method) {
          param_est <- M_step_variogram(
            eta_hat_list[[area_idx]],
            beta_tracker[, 1 + em_idx],
            TESSERAData_obj$X_list[[area_idx]],
            TESSERAData_obj$coords_list[[area_idx]],
            cov_type
          )
        } else if ("BRISC" == cov_fit_method) {
          param_est <- M_step_BRISC(
            eta_hat_list[[area_idx]],
            beta_tracker[, 1 + em_idx],
            TESSERAData_obj$X_list[[area_idx]],
            TESSERAData_obj$coords_list[[area_idx]],
            cov_type,
            nngp_k
          )
        } else {
          stop("Invalid cov_fit_method: 'BRISC' or 'variogram'")
        }
        
        
        cov_param_tracker[area_idx, 1 + em_idx, 1:length(param_est)] <- param_est
        
        # Update precision matrix Q and associated quantities
        param_est <- nngp_prec_mat(sp_dist_list[[area_idx]],
                                   TESSERAData_obj$coords_list[[area_idx]],
                                   cov_type,
                                   cov_param_tracker[area_idx, 1 + em_idx, ])
        Q_hat_list[[area_idx]] <- param_est$Q
        Dinv_list[[area_idx]] <- param_est$Dinv
        A_hat_list[[area_idx]] <- param_est$A
      }
      # Optimize in beta
      beta_tracker[, 1 + em_idx] <- M_step_beta(
        eta_hat_list,
        Q_hat_list,
        rep(1, n_areas),
        # tau^2
        TESSERAData_obj$X_list
      )[, 1]
    }
    
    if (verbose || (0 == (em_idx %% 100))) {
      message("beta ", paste(beta_tracker[, 1 + em_idx], collapse = " "), "\n")
      message("covariance parameters ",
              paste(cov_param_tracker[, 1 + em_idx, ], collapse = " "),
              "\n")
      message("log-lik ",
              paste(data_log_like_tracker[, em_idx], collapse = " "),
              "\n")
      
      message(
        paste(
          "Wrapping up of EM Iteration ",
          em_idx,
          " of ",
          em_iters,
          "; ",
          Sys.time() - t0_EM,
          " Elapsed"
        ),
        "\n"
      )
    }
    
    # Early stopping:
    #    NULL: Don't stop early.
    #    "abs_loglike": Difference in absolute data log likelihood, |old - new|.
    #    "rel_loglike": Relative difference in data log likelihood, |old - new|/|old|.
    #      Data log likelihood taken across alll areas.
    #    "abs_beta_norm": Absolute L2 norm of difference in beta, ||old - new||_2.
    #    "rel_beta_norm": Relative L2 norm of difference in beta,
    #      ||old - new||_2 / ||old||_2.
    if ((em_min_iters < em_idx) & !is.null(em_stopping)) {
      ll_old <- sum(data_log_like_tracker[, em_idx - 1], na.rm = TRUE)
      ll_new <- sum(data_log_like_tracker[, em_idx], na.rm = TRUE)
      beta_diff_norm <- sqrt(sum((beta_tracker[, em_idx + 1] - beta_tracker[, em_idx])^2))
      beta_old_norm <- sqrt(sum(beta_tracker[, em_idx]^2))
      if ("abs_loglike" == em_stopping) {
        if (abs(ll_new - ll_old) < em_tol) {
          message("Ending early", "\n")
          break
        }
      }
      else if ("rel_loglike" == em_stopping) {
        if (abs(ll_new - ll_old) / abs(ll_old) < em_tol) {
          message("Ending early", "\n")
          break
        }
      }
      else if ("abs_beta_norm" == em_stopping) {
        if (beta_diff_norm < em_tol) {
          message("Ending early", "\n")
          break
        }
      }
      else if ("rel_beta_norm" == em_stopping) {
        if (beta_diff_norm / beta_old_norm < em_tol) {
          message("Ending early", "\n")
          break
        }
      }
      else {
        warning("Invalid value for em_stopping.")
      }
    }
  }
  t1_EM <- Sys.time()
  message("Time ", t1_EM - t0_EM, "\n")
  
  # Compute Negative Hessians
  beta_neghessian <- neg_hessian_beta(Q_hat_list, rep(1, n_areas), TESSERAData_obj$X_list) # tau^2 is 1
  
  ## Name stuff
  
  rownames(beta_tracker) <- colnames(TESSERAData_obj$X_list[[1]])
  dimnames(cov_param_tracker) <- list(
    names(TESSERAData_obj$X_list),
    1:dim(cov_param_tracker)[2],
    c("Nugget", "Sill", "Range", "Smoothness")
  )
  rownames(R2_tracker) <- names(TESSERAData_obj$X_list)
  rownames(MSE_tracker) <- names(TESSERAData_obj$X_list)
  rownames(data_log_like_tracker) <- names(TESSERAData_obj$X_list)
  rownames(expected_log_like_tracker) <- names(TESSERAData_obj$X_list)
  dimnames(resid_moran) <- list(
    names(TESSERAData_obj$X_list),
    1:dim(resid_moran)[2],
    c("Moran_I", "ExpectedMoran_I", "PValue")
  )
  
  names(start_idx_list) <- names(TESSERAData_obj$X_list)
  
  rownames(beta_neghessian) <- colnames(TESSERAData_obj$X_list[[1]])
  colnames(beta_neghessian) <- colnames(TESSERAData_obj$X_list[[1]])
  
  rownames(fit_tracker) <- Reduce(c, lapply(TESSERAData_obj$counts_list, colnames))
  rownames(eta_tracker) <- rownames(fit_tracker)
  rownames(theta_tracker) <- rownames(fit_tracker)
  
  out <- (structure(
    list(
      # Coefficients, spatial parameters
      beta_hat = beta_tracker[, (em_idx + 1)],
      cov_param_hat = matrix(cov_param_tracker[, (em_idx + 1), ], nrow = n_areas, ncol = 4),
      
      # Fitted values and residuals, fitted Poisson parameters, estimated random effects
      predictions = fit_tracker[, em_idx],
      residuals = Reduce(c, z_list) - fit_tracker[, em_idx],
      eta_hat = eta_tracker[, em_idx],
      theta_hat = theta_tracker[, em_idx],
      phi_hat = eta_tracker[, em_idx] - Reduce(rbind, TESSERAData_obj$X_list) %*% beta_tracker[, (em_idx + 1)],
      
      # Full paths of coefficients, spatial parameters
      cov_param_tracker = cov_param_tracker[, 1:(em_idx + 1), ],
      beta_tracker = beta_tracker[, 1:(em_idx + 1)],
      
      # Trackers of various fit diagnostics
      R2_tracker = R2_tracker[, 1:em_idx],
      MSE_tracker = MSE_tracker[, 1:em_idx],
      data_log_like_tracker = data_log_like_tracker[, 1:em_idx],
      expected_log_like_tracker = expected_log_like_tracker[, 1:em_idx],
      resid_moran = resid_moran[, 1:em_idx, ],
      
      # Other utilities
      start_idx_list = start_idx_list,
      time = difftime(t1_EM, t0_EM),
      
      # Hessians (for standard errors)
      beta_neghessian = beta_neghessian,
      
      run_settings = list(
        gene_name = as.character(gene_name),
        gene_idx = gene_idx,
        model_type = "spNNGP",
        cov_type = as.character(cov_type),
        nngp_k = nngp_k,
        em_iters = em_iters,
        opt_iters = opt_iters,
        em_min_iters = em_min_iters,
        em_tol = em_tol,
        em_stopping = em_stopping,
        beta_init = beta_init,
        cov_init = cov_init,
        cov_fit_method = cov_fit_method,
        verbose = verbose,
        em_iters_actual = em_idx
      )
    ),
    class = "TESSERAOutput"
  ))
  
  out$performanceSummary <- summarizeTESSERAPerformance(TESSERAData_obj, out)
  return(out)
}

#' Check inputs for the TESSERA_spNNGP method.
#'  Make sure that the input object has everything needed to run without error.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param TESSERAData_obj Object containing data.
#'  Created by the prepData method.
#'
#' @note Does not return anything.
#' @note This method can be used to check a hand-created input object.
#'  E.g., if a user does not want to use prepData.
#'
#' @returns Nothing.
#' @import Matrix
#' 
#' @export
#' 
#' @example inst/examples/gen_data_and_run.R
checkInputsTESSERAspNNGP <- function (TESSERAData_obj) {
  # Check that the bare minimum is present
  stopifnot(!is.null(TESSERAData_obj$counts_list))
  stopifnot(!is.null(TESSERAData_obj$X_list))
  stopifnot(!is.null(TESSERAData_obj$coords_list))
  stopifnot(!is.null(TESSERAData_obj$library_size_list))
  
  # Counts-specific checks
  # Check that there is at least one gene, really that it's a matrix.
  stopifnot(1 <= min(sapply(TESSERAData_obj$counts_list, nrow)))
  # Check that same number of genes present
  stopifnot(1 == length(unique(
    lapply(TESSERAData_obj$counts_list, nrow)
  )))
  # Check ordering by rownames (gene names)
  stopifnot(all(sapply(
    lapply(TESSERAData_obj$counts_list, rownames),
    identical,
    rownames(TESSERAData_obj$counts_list[[1]])
  )))
  
  # Coordinates-specific checks
  # Need at least 2 coordinates
  stopifnot(1 < min(sapply(TESSERAData_obj$coords_list, ncol)))
  
  # Check that the number of measurements/cells is the same across lists
  # Implicit check that the number of entries in list is the same (number of samples)
  stopifnot(identical(
    sapply(TESSERAData_obj$counts_list, ncol),
    sapply(TESSERAData_obj$X_list, nrow)
  ))
  stopifnot(identical(
    sapply(TESSERAData_obj$counts_list, ncol),
    sapply(TESSERAData_obj$library_size_list, length)
  ))
  stopifnot(identical(
    sapply(TESSERAData_obj$counts_list, ncol),
    sapply(TESSERAData_obj$coords_list, nrow)
  ))
  
  # Check that the names of measurements match
  for (idx in 1:length(TESSERAData_obj$counts_list)) {
    stopifnot(identical(
      colnames(TESSERAData_obj$counts_list[[idx]]),
      rownames(TESSERAData_obj$X_list[[idx]])
    ))
    stopifnot(identical(
      colnames(TESSERAData_obj$counts_list[[idx]]),
      names(TESSERAData_obj$library_size_list[[idx]])
    ))
    stopifnot(identical(
      colnames(TESSERAData_obj$counts_list[[idx]]),
      rownames(TESSERAData_obj$coords_list[[idx]])
    ))
  }
  
  # Check ordering of lists (sample IDs)
  stopifnot(identical(
    names(TESSERAData_obj$counts_list),
    names(TESSERAData_obj$X_list)
  ))
  stopifnot(identical(
    names(TESSERAData_obj$counts_list),
    names(TESSERAData_obj$library_size_list)
  ))
  stopifnot(identical(
    names(TESSERAData_obj$counts_list),
    names(TESSERAData_obj$coords_list)
  ))
}
