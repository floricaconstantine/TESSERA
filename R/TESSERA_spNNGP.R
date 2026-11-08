## Wrapper functions to fit model.
# Dependencies in file: Matrix.
# Dependencies: Functions from utils.R, models.R, E_step.R, M_step.R.
# Dependences from functions in other files not listed: pracma, gstat, BRISC.
# Rcpp dependencies: calc_moran.cpp, nngp_prec_mat.cpp.


#' Fit Multi-Sample Poisson Spatial GLMM via spNNGP
#'
#' Fits a multi-sample Poisson spatial generalized linear mixed model (GLMM)
#' using a shared set of fixed effects across all samples while permitting
#' sample-specific spatial random effects modeled via a Sparse Nearest Neighbor
#' Gaussian Process (spNNGP). Supports Exponential, Matern, Gaussian, and
#' Spherical covariance kernels.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param TESSERAData_obj An object containing prepared data, typically
#'   created by \code{\link{prep_data}}.
#' @param gene_name Character: The name of the gene/measurement (row) to fit.
#' @param cov_type Character: The covariance kernel for the Gaussian process.
#'   Options are "Exp", "Mat", "Gau", and "Sph".
#' @param nngp_k Integer: The number of nearest neighbors to use for the
#'   spNNGP kernel approximation.
#' @param em_iters Integer: Maximum number of ECM iterations.
#' @param opt_iters Integer: Number of inner CM steps (Conditional Maximization)
#'   per iteration. The algorithm maximizes the expected likelihood for
#'   the covariance parameters, then \eqn{\beta}, holding other parameters constant.
#' @param em_min_iters Integer: Minimum number of ECM iterations to perform
#'   before allowing early stopping.
#' @param em_tol Numeric: Convergence tolerance for early stopping.
#' @param em_stopping Character: Metric used for early stopping:
#' \itemize{
#'   \item \code{NULL}: No early stopping.
#'   \item "abs_loglike": Absolute change in total log-likelihood.
#'   \item "rel_loglike": Relative change in total log-likelihood.
#'   \item "abs_beta_norm": Absolute \eqn{L_2} norm of the change in \eqn{\beta}.
#'   \item "rel_beta_norm": Relative \eqn{L_2} norm of the change in \eqn{\beta}.
#' }
#' @param beta_init Initial value for \eqn{\beta}. Options: "glm" (fit a Poisson GLM),
#'   "random" (standard normal), or a numeric vector.
#' @param cov_init Method for initializing covariance parameters. Options:
#'   "BRISC" (calls BRISC on residuals), "variogram" (initializes via empirical
#'   variogram), or a numeric vector/matrix.
#' @param cov_fit_method Character: The method used for the M-step update
#'   of covariance parameters. Options are "BRISC" or "variogram".
#' @param verbose Logical: Whether to print iteration-wise parameter updates.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \strong{beta_hat}: Estimated fixed effect coefficients.
#'   \item \strong{cov_param_hat}: Estimated spatial covariance parameters
#'     (nugget, sill, range, and optionally smoothness).
#'   \item \strong{phi_hat}: Estimated spatial random effects.
#'   \item \strong{theta_hat}: Estimated Poisson rate parameters \eqn{exp(X\beta + \phi)}.
#'   \item \strong{beta_tracker}: Matrix of \eqn{\beta} estimates across iterations.
#'   \item \strong{cov_param_tracker}: Array of covariance parameter history
#'     (areas x iterations x parameters).
#'   \item \strong{data_log_like_tracker}: Total log-likelihood across iterations.
#'   \item \strong{MSE_tracker}: Mean Squared Error history across iterations.
#'   \item \strong{beta_neghessian}: Negative Hessian of the final \eqn{\beta} estimate.
#'   \item \strong{performanceSummary}: A summary data frame for each sample;
#'     see \code{\link{summarize_TESSERA}}.
#'   \item \strong{time}: Total execution time.
#'   \item \strong{run_settings}: List of parameters used for the run.
#' }
#'
#' @references Meng, Xiao-Li, and Donald B. Rubin. "Maximum likelihood estimation via the ECM algorithm: A general framework." Biometrika 80.2 (1993): 267-278.
#' @references Saha, Arkajyoti, and Abhirup Datta. "BRISC: bootstrap for rapid inference on spatial covariances." Stat 7.1 (2018): e184.
#'
#' @note The spNNGP approach is designed for scalability in datasets where
#'   traditional lattice adjacency is difficult to define.
#'
#' @import Matrix
#' @importFrom Matrix diag
#' @importFrom sparseinv cholPermute
#' @importFrom stats coef cor dpois poisson predict rnorm var
#' @importFrom Rcpp sourceCpp evalCpp
#' @useDynLib TESSERA
#' @export
#'
#' @examples
#' # Locate the prepped TESSERAData object in inst/extdata
#' rds_path <- system.file("extdata", "example_prepData.rds", package = "TESSERA")
#' # Load the TESSERAData object
#' TESSERA_data <- readRDS(rds_path)
#'
#' # Fit the Poisson generalized spatial linear model using spNNGP
#' TESSERA_out_spNNGP <- suppressMessages(suppressWarnings(
#'   TESSERA_spNNGP(
#'     TESSERAData_obj = TESSERA_data,
#'     gene_name = "example",
#'     cov_type = "Exp",
#'     em_iters = 2,
#'     opt_iters = 1,
#'     verbose = FALSE
#'   )
#' ))
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
                           verbose = FALSE)
{
  # Start clock
  t0_EM <- Sys.time()
  
  # Check inputs
  check_inputs_TESSERA_spNNGP(TESSERAData_obj)
  
  # Extract gene of interest and associated counts
  gene_idx <- which(rownames(TESSERAData_obj$counts_list[[1]]) == gene_name)
  z_list <- lapply(TESSERAData_obj$counts_list, function (x) {
    x[gene_idx, ]
  })
  
  # Number of areas
  n_areas <- length(z_list)
  
  # ==========================================
  # --- Spatial Coordinate Sorting Block ---
  perm_list <- list()
  rev_perm_list <- list()
  
  for (area_idx in 1:n_areas) {
    coords <- TESSERAData_obj$coords_list[[area_idx]]
    
    # Create the sorting permutation (sweep left-to-right, bottom-to-top)
    perm <- order(coords[, 1], coords[, 2])
    
    # Save the permutation and the REVERSE permutation
    perm_list[[area_idx]] <- perm
    rev_perm_list[[area_idx]] <- order(perm)
    
    # Apply permutation to all relevant data
    TESSERAData_obj$coords_list[[area_idx]] <- coords[perm, , drop = FALSE]
    TESSERAData_obj$X_list[[area_idx]] <- TESSERAData_obj$X_list[[area_idx]][perm, , drop = FALSE]
    TESSERAData_obj$library_size_list[[area_idx]] <- TESSERAData_obj$library_size_list[[area_idx]][perm]
    z_list[[area_idx]] <- z_list[[area_idx]][perm]
  }
  # ---------------------------------------------
  # ==========================================
  
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
    sp_dist_list[[idx]] <- sparse_dist(TESSERAData_obj$coords_list[[idx]], nngp_k)
  }
  
  # ==========================================
  # --- Precompute Static Structures Outside the EM Loop ---
  # 1. Static Neighbor Distance Matrices (in C++)
  # 2. Moran's I Coordinate Constants (W, S1, S2, ei)
  # 3. E-step etahat static term1 and 0.5 + z
  nb_dist_list <- vector("list", n_areas)
  moran_precomp_list <- vector("list", n_areas)
  etahat_term1_list <- vector("list", n_areas)
  z_half_list <- vector("list", n_areas)
  
  for (area_idx in 1:n_areas) {
    coords_area <- as.matrix(TESSERAData_obj$coords_list[[area_idx]])
    
    # Fast C++ neighbor distance matrix precomputation
    nb_dist_list[[area_idx]] <- build_nb_dist_cpp(coords_area, sp_dist_list[[area_idx]], nngp_k)
    
    # Precompute coordinate-only spatial constants for Moran's I
    moran_precomp_list[[area_idx]] <- precomp_moran(coords_area[, 1], coords_area[, 2])
    
    # Precompute 0.5 + z and static term1 for E_step_etahat
    t1_tmp <- as.numeric(z_list[[area_idx]] + 0.5)
    z_half_list[[area_idx]] <- t1_tmp
    etahat_term1_list[[area_idx]] <- t1_tmp * log(t1_tmp / TESSERAData_obj$library_size_list[[area_idx]]) - 0.5
  }
  # ---------------------------------------------
  # ==========================================
  
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
    # Stack into a single vector/matrix efficiently
    z_vec <- unlist(z_list, use.names = FALSE)
    lib_vec <- unlist(TESSERAData_obj$library_size_list, use.names = FALSE)
    X_mat <- do.call(rbind, TESSERAData_obj$X_list)
    
    # Fit GLM
    beta_tmp <- as.vector(stats::coef(
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
                               nb_dist_list[[area_idx]],
                               cov_type,
                               cov_param_tracker[area_idx, 1, ])
    Q_hat_list[[area_idx]] <- param_est$Q
    Dinv_list[[area_idx]] <- param_est$Dinv
    A_hat_list[[area_idx]] <- param_est$A
  }
  
  # Store permutation pattern
  V_chol_perm_list <- list()
  
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
    eta_hat_list <- list()
    
    for (area_idx in 1:n_areas) {
      Vinv <- build_Vinv_cpp(Q_hat_list[[area_idx]], 1.0, z_half_list[[area_idx]])
      
      # Compute Cholesky factorization once per E-step (and cache permutation P on iter 1)
      if (1 == em_idx) {
        chol_current <- sparseinv::cholPermute(Vinv)
        V_chol_perm_list[[area_idx]] <- list(P = chol_current$P)
      } else {
        P_mat <- V_chol_perm_list[[area_idx]]$P
        chol_current <- list(Qpermchol = Matrix::t(Matrix::chol(
          Matrix::forceSymmetric(Matrix::t(P_mat) %*% Vinv %*% P_mat)
        )),
        P = P_mat)
      }
      
      # Compute Vhat reusing chol_current$Qpermchol (avoids 2nd cholPermute inside Takahashi_Davis)
      Vhat_current <- E_step_Vhat(
        Vinv,
        P = chol_current$P,
        cholQp = chol_current$Qpermchol,
        return_dgC = TRUE
      )
      
      # Get the mean reusing chol_current (solves via triangular L, eliminating sparse LU)
      eta_hat_list[[area_idx]] <- E_step_etahat(
        Vinv,
        Q_hat_list[[area_idx]],
        1.0,
        beta_tracker[, em_idx],
        TESSERAData_obj$X_list[[area_idx]],
        z_list[[area_idx]],
        TESSERAData_obj$library_size_list[[area_idx]],
        term1_precomp = etahat_term1_list[[area_idx]],
        chol_perm = chol_current
      )
      rm(chol_current)
      
      # Get a few more things out of the E-step: theta and predictions
      theta_hat <- as.numeric(E_step_thetahat(Vhat_current, eta_hat_list[[area_idx]]))
      z_hat <- as.numeric(E_step_predict(theta_hat, TESSERAData_obj$library_size_list[[area_idx]]))
      
      # Store stuff
      fit_tracker[start_idx_list[area_idx]:(start_idx_list[area_idx] + length(z_hat) - 1), em_idx] <- z_hat
      eta_tracker[start_idx_list[area_idx]:(start_idx_list[area_idx] + length(eta_hat_list[[area_idx]]) - 1), em_idx] <- as.numeric(eta_hat_list[[area_idx]])
      theta_tracker[start_idx_list[area_idx]:(start_idx_list[area_idx] + length(theta_hat) - 1), em_idx] <- theta_hat
      
      # Performance
      R2_tracker[area_idx, em_idx] <- stats::cor(z_hat, z_list[[area_idx]])^2
      MSE_tracker[area_idx, em_idx] <- mean(abs(z_hat - z_list[[area_idx]])^2)
      
      # Observed, expected, sd (using precomputed coordinate constants)
      resid_moran[area_idx, em_idx, ] <- calc_moran_precomp(
        z_hat - z_list[[area_idx]],
        TESSERAData_obj$coords_list[[area_idx]][, 1],
        TESSERAData_obj$coords_list[[area_idx]][, 2],
        moran_precomp_list[[area_idx]]
      )
      
      # Data log likelihood
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
        Vhat_current,
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
      
      # DESTROY the sparse Vhat to prevent Memory Overflow
      rm(Vhat_current)
    }
    
    # Run (C)M-Step: Optimize
    for (opt_idx in 1:opt_iters) {
      # Initialize with current values (needed)
      if (1 == opt_idx) {
        beta_tracker[, 1 + em_idx] <- beta_tracker[, em_idx]
      }
      
      # Only build the sparse matrix Q on the final inner CM step (opt_idx == opt_iters);
      # on earlier inner steps, compute X^T Q eta and X^T Q X directly in C++.
      is_last_opt <- (opt_idx == opt_iters)
      XtQX_step_list <- if (!is_last_opt)
        vector("list", n_areas)
      else
        NULL
      XtQeta_step_list <- if (!is_last_opt)
        vector("list", n_areas)
      else
        NULL
      
      for (area_idx in 1:n_areas) {
        # Initialize with current values (just as a placeholder)
        if (1 == opt_idx) {
          cov_param_tracker[area_idx, 1 + em_idx, ] <- cov_param_tracker[area_idx, em_idx, ]
        }
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
            nngp_k,
            prev_params = NULL # cov_param_tracker[area_idx, 1 + em_idx, ]
          )
        } else {
          stop("Invalid cov_fit_method: 'BRISC' or 'variogram'")
        }
        
        cov_param_tracker[area_idx, 1 + em_idx, 1:length(param_est)] <- param_est
        
        # Update precision matrix Q (on final opt_idx) or direct beta cross-products (on opt_idx < opt_iters)
        prec_out <- nngp_prec_mat(
          sp_dist_list[[area_idx]],
          nb_dist_list[[area_idx]],
          cov_type,
          cov_param_tracker[area_idx, 1 + em_idx, ],
          build_Q = is_last_opt,
          X = TESSERAData_obj$X_list[[area_idx]],
          eta = eta_hat_list[[area_idx]]
        )
        if (is_last_opt) {
          Q_hat_list[[area_idx]] <- prec_out$Q
          Dinv_list[[area_idx]] <- prec_out$Dinv
          A_hat_list[[area_idx]] <- prec_out$A
        } else {
          Dinv_list[[area_idx]] <- prec_out$Dinv
          XtQX_step_list[[area_idx]] <- prec_out$XtQX
          XtQeta_step_list[[area_idx]] <- prec_out$XtQeta
        }
      }
      # Optimize in beta
      beta_tracker[, 1 + em_idx] <- M_step_beta(
        eta_list = eta_hat_list,
        Q_list = Q_hat_list,
        tau2_list = rep(1, n_areas),
        # tau^2
        X_list = TESSERAData_obj$X_list,
        model_type = "spNNGP",
        XtQX_list = XtQX_step_list,
        XtQeta_list = XtQeta_step_list
      )
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
  
  # Compute Negative Hessians (MUST HAPPEN BEFORE UN-SORTING)
  beta_neghessian <- neg_hessian_beta(Q_hat_list, rep(1, n_areas), TESSERAData_obj$X_list) # tau^2 is 1
  
  # ==========================================
  # --- Reverse Permutation Block ---
  # Un-sort the trackers and input data so output matches the original counts_list
  for (area_idx in 1:n_areas) {
    # Get the row indices for this specific area in the flat trackers
    idx_range <- start_idx_list[area_idx]:(start_idx_list[area_idx] + length(z_list[[area_idx]]) - 1)
    
    # Apply the reverse order
    rev_perm <- rev_perm_list[[area_idx]]
    
    # Un-sort EM trackers
    fit_tracker[idx_range, ] <- fit_tracker[idx_range, , drop = FALSE][rev_perm, , drop = FALSE]
    eta_tracker[idx_range, ] <- eta_tracker[idx_range, , drop = FALSE][rev_perm, , drop = FALSE]
    theta_tracker[idx_range, ] <- theta_tracker[idx_range, , drop = FALSE][rev_perm, , drop = FALSE]
    
    # Un-sort the data lists so formulas in the output block align perfectly
    z_list[[area_idx]] <- z_list[[area_idx]][rev_perm]
    TESSERAData_obj$X_list[[area_idx]] <- TESSERAData_obj$X_list[[area_idx]][rev_perm, , drop = FALSE]
    TESSERAData_obj$coords_list[[area_idx]] <- TESSERAData_obj$coords_list[[area_idx]][rev_perm, , drop = FALSE]
    TESSERAData_obj$library_size_list[[area_idx]] <- TESSERAData_obj$library_size_list[[area_idx]][rev_perm]
  }
  # ----------------------------------------
  # ==========================================
  
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
  
  rownames(fit_tracker) <- unlist(lapply(TESSERAData_obj$counts_list, colnames), use.names = FALSE)
  rownames(eta_tracker) <- rownames(fit_tracker)
  rownames(theta_tracker) <- rownames(fit_tracker)
  
  out <- (structure(
    list(
      # Coefficients, spatial parameters
      beta_hat = beta_tracker[, (em_idx + 1)],
      cov_param_hat = matrix(cov_param_tracker[, (em_idx + 1), ], nrow = n_areas, ncol = 4),
      
      # Fitted values and residuals, fitted Poisson parameters, estimated random effects
      predictions = fit_tracker[, em_idx],
      residuals = unlist(z_list, use.names = FALSE) - fit_tracker[, em_idx],
      eta_hat = eta_tracker[, em_idx],
      theta_hat = theta_tracker[, em_idx],
      phi_hat = eta_tracker[, em_idx] - do.call(rbind, TESSERAData_obj$X_list) %*% beta_tracker[, (em_idx + 1)],
      
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
  
  out$performanceSummary <- summarize_TESSERA(TESSERAData_obj, out)
  return(out)
}



#' Check inputs for the TESSERA_spNNGP method.
#'  Make sure that the input object has everything needed to run without error.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param TESSERAData_obj Object containing data.
#'  Created by the prep_data method.
#'
#' @note Does not return anything.
#' @note This method can be used to check a hand-created input object.
#'  E.g., if a user does not want to use prep_data.
#'
#' @returns Nothing.
#' @import Matrix
check_inputs_TESSERA_spNNGP <- function (TESSERAData_obj) {
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
