## M-Step (Maximization Step) functions.
# Dependencies in file: Matrix, pracma, stats, gstat, BRISC.

#' Maximize the expected likelihood in tau^2, holding other variables constant.
#' Part of the M-Step in the EM algorithm.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param trace_scalars List of precomputed trace scalars from the E-step.
#' @param gamma_val Current estimate of correlation parameter.
#' @param model_type String for the spatial model ("CAR", "SAR", or "Leroux").
#' @param eta_hat Estimated mean of eta.
#' @param Q Unscaled precision matrix.
#' @param beta_hat Current estimate of covariate effects.
#' @param X Covariate matrix.
#'
#' @returns Maximum likelihood estimate of scaling parameter tau^2.
M_step_tau2 <- function(trace_scalars,
                        gamma_val,
                        model_type,
                        eta_hat,
                        Q,
                        beta_hat,
                        X) {
  # eta - X beta
  vector_term <- eta_hat - (X %*% beta_hat)
  
  # (eta - X beta)^\top Q (eta - X beta)
  # Speed: crossprod is optimized in R
  term1 <- as.numeric(crossprod(Q %*% vector_term, vector_term))
  
  # Trace[Q V] calculated via trace decomposition
  # Reconstructs the exact trace instantly using O(1) memory
  if ("CAR" == model_type) {
    # Q = D - gamma * W
    term2 <- trace_scalars$tr_DV - (gamma_val * trace_scalars$tr_WV)
  } else if ("SAR" == model_type) {
    # Q = D - 2*gamma*W + gamma^2 * W Z
    term2 <- trace_scalars$tr_DV - (2.0 * gamma_val * trace_scalars$tr_WV) + ((gamma_val^2) * trace_scalars$tr_WZV)
  } else if ("Leroux" == model_type) {
    # Q = I + gamma * (D - W - I)
    term2 <- trace_scalars$tr_V + (gamma_val * trace_scalars$tr_DWIV)
  } else {
    stop("Invalid model_type in M_step_tau2")
  }
  
  # tau^2 = (1/n) (term1 + term2)
  tau2_hat <- (term1 + term2) / dim(Q)[1]
  
  # Make sure tau^2 doesn't go off the rails
  close_to_zero_const <- 1e2 * .Machine$double.eps
  if (close_to_zero_const >= tau2_hat) {
    warning("Invalid tau^2, setting to a small positive number")
    tau2_hat <- max(tau2_hat, close_to_zero_const)
  }
  
  return(tau2_hat)
}


#' M-step optimization for beta (covariate coefficients).
#' Part of the M-Step in the EM algorithm.
#'
#' @param eta_list List of estimated means of eta.
#' @param Q_list List of unscaled precision matrices.
#' @param tau2_list List of precision/covariance matrix scalings.
#' @param X_list List of covariate matrices.
#' @param model_type "CAR", "SAR", "Leroux", or "spNNGP".
#' @param gamma_list List of current gamma estimates (required for lattice).
#' @param XtDX_list Precomputed X^T D X (required for CAR/SAR).
#' @param XtWX_list Precomputed X^T W X (required for CAR/SAR).
#' @param XtWZX_list Precomputed X^T W Z X (required for SAR).
#' @param XtX_list Precomputed X^T X (required for Leroux).
#' @param XtDWIX_list Precomputed X^T (D - W - I) X (required for Leroux).
#'
#' @returns Estimated beta vector.
#'
#' @importFrom Matrix solve crossprod
M_step_beta <- function(eta_list, Q_list, tau2_list, X_list,
                        model_type = "spNNGP",
                        gamma_list = NULL,
                        XtDX_list = NULL,
                        XtWX_list = NULL,
                        XtWZX_list = NULL,
                        XtX_list = NULL,
                        XtDWIX_list = NULL) {
  
  n_areas <- length(eta_list)
  p <- ncol(X_list[[1]])
  
  B <- matrix(0, nrow = p, ncol = p)
  zeta_vec <- matrix(0, nrow = p, ncol = 1)
  
  for (idx in 1:n_areas) {
    # 1. Vector part: X^T Q eta 
    # Dynamic, but incredibly fast O(Np) operation
    Q_inv_tau2 <- Q_list[[idx]] / tau2_list[[idx]]
    XtQ <- Matrix::crossprod(X_list[[idx]], Q_inv_tau2)
    zeta_vec <- zeta_vec + as.numeric(XtQ %*% eta_list[[idx]])
    
    # 2. Matrix part: X^T Q X 
    # Assembled from precomputed matrices for Lattice, dynamic for spNNGP
    if (model_type == "CAR") {
      B_area <- (XtDX_list[[idx]] - gamma_list[idx] * XtWX_list[[idx]]) / tau2_list[[idx]]
      B <- B + B_area
      
    } else if (model_type == "SAR") {
      B_area <- (XtDX_list[[idx]] - 2 * gamma_list[idx] * XtWX_list[[idx]] + 
                   (gamma_list[idx]^2) * XtWZX_list[[idx]]) / tau2_list[[idx]]
      B <- B + B_area
      
    } else if (model_type == "Leroux") {
      B_area <- (gamma_list[idx] * XtDWIX_list[[idx]] + XtX_list[[idx]]) / tau2_list[[idx]]
      B <- B + B_area
      
    } else if (model_type == "spNNGP") {
      # Dynamic fallback for spatial covariance matrices
      B <- B + as.matrix(XtQ %*% X_list[[idx]])
      
    } else {
      stop("Invalid model type in M_step_beta")
    }
  }
  
  # beta = B^{-1} zeta
  # Force back to base R to avoid S4 dispatch errors, beta = B^{-1} zeta
  beta_hat <- base::solve(as.matrix(B), as.numeric(zeta_vec))
  return(as.numeric(beta_hat))
}


## Model-specific updates for gamma

#' Maximize the expected likelihood in gamma, holding other variables constant.
#' Part of the M-Step in the EM algorithm.
#' ONLY APPLIES TO THE CAR MODEL.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param trace_scalars List of precomputed trace scalars from the E-step.
#' @param eta_hat Estimated mean of eta.
#' @param tau2 Precision matrix scaling parameter tau^2.
#' @param beta_hat Current estimate of covariate effects.
#' @param X Covariate matrix.
#' @param W Neighbor/adjacency matrix (symmetric, binary).
#' @param D Diagonal degree matrix (row sums of W).
#'  Not needed, but passed in for consistency with API.
#' @param eig_vals Eigenvalues of Z = D^\{-1\} W.
#' @param gamma_current Previous value of gamma.
#'  Only needed for failure modes.
#'
#' @returns Maximum likelihood estimate of correlation parameter gamma.
#'  A list with gamma_hat (estimate) and grad_val (gradient).
#'
#' @note Requires the pracma library.
#'
#' @importFrom pracma newtonRaphson
#' @importFrom pracma fzero
#' @importFrom stats uniroot
#' @importFrom stats runif
M_step_gamma_CAR <- function(trace_scalars,
                             eta_hat,
                             tau2,
                             beta_hat,
                             X,
                             W,
                             D,
                             eig_vals,
                             gamma_current = NULL) {
  # eta - X beta
  vector_term <- eta_hat - (X %*% beta_hat)
  
  # (eta - X beta)^\top W (eta - X beta)
  term1 <- as.numeric(crossprod(W %*% vector_term, vector_term))
  
  # Tr[W V]
  term2 <- trace_scalars$tr_WV
  
  # (term1 + term2) / (2 tau^2)
  constant_term <- (term1 + term2) * (0.5 / tau2)
  
  # Gradient function
  grad_fcn <- function(x) {
    # (-1/2) sum_i lambda_i / (1 - gamma lambda_i)
    eig_term <- (-0.5) * sum(eig_vals / (1.0 - (x * eig_vals)))
    return (eig_term + constant_term)
  }
  
  # Set upper bound
  if (0 == sum(abs(eig_vals - 1) < 1e-4)) {
    upper <- 1.0 # No eigenvalues of 1.0
  } else {
    upper <- 1.0 - 1e-3
  }
  
  # Try basic uniroot first
  err_flag <- tryCatch({
    # Find root
    gamma_out <- stats::uniroot(
      grad_fcn,
      lower = 0.0,
      # -1.0,
      upper = upper,
      extendInt = "yes"
    )
    
    # Return estimate and gradient value
    return(list(gamma_hat = gamma_out$root, grad_val = gamma_out$f.root))
  }, error = function(cond) {
    warning("UNIROOT FAILED")
    err_flag <- TRUE
  })
  # Then try Newton-Raphson
  if (err_flag) {
    # If null, random
    if (is.null(gamma_current)) {
      gamma_current <- stats::runif(1)
      message("Assigning random gamma ", gamma_current, "\n")
    }
    
    err_flag <- tryCatch({
      gamma_out <- pracma::newtonRaphson(grad_fcn, gamma_current)
      
      # Return estimate and gradient value
      return(list(gamma_hat = gamma_out$root, grad_val = gamma_out$f.root))
    }, error = function(cond) {
      warning("NEWTON-RAPHSON FAILED")
      err_flag <- TRUE
    })
  }
  # Then try fzero
  if (err_flag) {
    # If null, random
    if (is.null(gamma_current)) {
      gamma_current <- stats::runif(1)
      message("Assigning random gamma ", gamma_current, "\n")
    }
    
    err_flag <- tryCatch({
      gamma_out <- pracma::fzero(grad_fcn, gamma_current)
      
      # Return estimate and gradient value
      return(list(gamma_hat = gamma_out$x, grad_val = gamma_out$fval))
    }, error = function(cond) {
      warning("FZERO FAILED")
      err_flag <- TRUE
    })
  }
  if (err_flag) {
    stop("OPTIMIZATION FOR GAMMA FAILED: CANNOT CONTINUE")
  }
}


#' Maximize the expected likelihood in gamma, holding other variables constant.
#' Part of the M-Step in the EM algorithm.
#' ONLY APPLIES TO THE SAR MODEL.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param trace_scalars List of precomputed trace scalars from the E-step.
#' @param eta_hat Estimated mean of eta.
#' @param tau2 Precision matrix scaling parameter tau^2.
#' @param beta_hat Current estimate of covariate effects.
#' @param X Covariate matrix.
#' @param W Neighbor/adjacency matrix (symmetric, binary).
#' @param D Diagonal degree matrix (row sums of W).
#' @param eig_vals Eigenvalues of Z = D^\{-1\} W.
#' @param gamma_current Previous value of gamma.
#'  Only needed for failure modes.
#'
#' @returns Maximum likelihood estimate of correlation parameter gamma.
#'  A list with gamma_hat (estimate) and grad_val (gradient).
#'
#' @note Requires the Matrix library.
#' @note Requires the pracma library.
#'
#' @import Matrix
#' @importFrom Matrix Diagonal
#' @importFrom Matrix diag
#' @importFrom pracma newtonRaphson
#' @importFrom pracma fzero
#' @importFrom stats uniroot
#' @importFrom stats runif
M_step_gamma_SAR <- function(trace_scalars,
                             eta_hat,
                             tau2,
                             beta_hat,
                             X,
                             W,
                             D,
                             eig_vals,
                             gamma_current = NULL) {
  # D^{-1}
  D_inv <- Matrix::Diagonal(dim(W)[1], 1 / Matrix::diag(D))
  # Z = D^{-1} W
  Z <- D_inv %*% W
  
  # eta - X beta
  vector_term <- eta_hat - (X %*% beta_hat)
  
  # (eta - X beta)^\top (2 W) (eta - X beta)
  zeta1 <- 2.0 * as.numeric(crossprod(W %*% vector_term, vector_term))
  
  # (eta - X beta)^\top W Z (eta - X beta)
  zeta2 <- as.numeric(crossprod(W %*% vector_term, Z %*% vector_term))
  
  # Tr[W V]
  zeta3 <- trace_scalars$tr_WV
  
  # Tr[W Z V]
  zeta4 <- trace_scalars$tr_WZV
  
  grad_fcn <- function(x) {
    # (-1) sum_i lambda_i / (1 - gamma lambda_i)
    eig_term <- (-1.0) * sum(eig_vals / (1.0 - (x * eig_vals)))
    
    # (-1/2 tau^2) (-zeta1 + 2 gamma zeta2)
    term1 <- (-0.5 / tau2) * (-zeta1 + x * (2.0 * zeta2))
    # (-1/tau^2) (-zeta3 + gamma zeta4)
    term2 <- (-1.0 / tau2) * (-zeta3 + x * zeta4)
    
    return(eig_term + term1 + term2)
  }
  
  # Set upper bound
  if (0 == sum(abs(eig_vals - 1) < 1e-4)) {
    upper <- 1.0 # No eigenvalues of 1.0
  } else {
    upper <- 1.0 - 1e-3
  }
  
  # Try basic uniroot first
  err_flag <- tryCatch({
    # Find root
    gamma_out <- stats::uniroot(
      grad_fcn,
      lower = 0.0,
      # -1.0,
      upper = upper,
      extendInt = "yes"
    )
    
    # Return estimate and gradient value
    return(list(gamma_hat = gamma_out$root, grad_val = gamma_out$f.root))
  }, error = function(cond) {
    warning("UNIROOT FAILED")
    err_flag <- TRUE
  })
  # Then try Newton-Raphson
  if (err_flag) {
    # If null, random
    if (is.null(gamma_current)) {
      gamma_current <- stats::runif(1)
      message("Assigning random gamma ", gamma_current, "\n")
    }
    
    err_flag <- tryCatch({
      gamma_out <- pracma::newtonRaphson(grad_fcn, gamma_current)
      
      # Return estimate and gradient value
      return(list(gamma_hat = gamma_out$root, grad_val = gamma_out$f.root))
    }, error = function(cond) {
      warning("NEWTON-RAPHSON FAILED")
      err_flag <- TRUE
    })
  }
  # Then try fzero
  if (err_flag) {
    # If null, random
    if (is.null(gamma_current)) {
      gamma_current <- stats::runif(1)
      message("Assigning random gamma ", gamma_current, "\n")
    }
    
    err_flag <- tryCatch({
      gamma_out <- pracma::fzero(grad_fcn, gamma_current)
      
      # Return estimate and gradient value
      return(list(gamma_hat = gamma_out$x, grad_val = gamma_out$fval))
    }, error = function(cond) {
      warning("FZERO FAILED")
      err_flag <- TRUE
    })
  }
  if (err_flag) {
    stop("OPTIMIZATION FOR GAMMA FAILED: CANNOT CONTINUE")
  }
}


#' Maximize the expected likelihood in gamma, holding other variables constant.
#' Part of the M-Step in the EM algorithm.
#' ONLY APPLIES TO THE LEROUX MODEL.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param trace_scalars List of precomputed trace scalars from the E-step.
#' @param eta_hat Estimated mean of eta.
#' @param tau2 Precision matrix scaling parameter tau^2.
#' @param beta_hat Current estimate of covariate effects.
#' @param X Covariate matrix.
#' @param W Neighbor/adjacency matrix (symmetric, binary).
#' @param D Diagonal degree matrix (row sums of W).
#' @param eig_vals Eigenvalues of Z = D - W.
#' @param gamma_current Previous value of gamma.
#'  Only needed for failure modes.
#'
#' @returns Maximum likelihood estimate of correlation parameter gamma.
#'  A list with gamma_hat (estimate) and grad_val (gradient).
#'
#' @note Requires the Matrix library.
#'
#' @import Matrix
#' @importFrom Matrix Diagonal
#' @importFrom pracma newtonRaphson
#' @importFrom pracma fzero
#' @importFrom stats uniroot
#' @importFrom stats runif
M_step_gamma_Leroux <- function(trace_scalars,
                                eta_hat,
                                tau2,
                                beta_hat,
                                X,
                                W,
                                D,
                                eig_vals,
                                gamma_current = NULL) {
  # Identity matrix
  id_mat <- Matrix::Diagonal(dim(W)[1], 1)
  # D - W - I
  DWI <- D - W - id_mat
  
  # eta - X beta
  vector_term <- eta_hat - (X %*% beta_hat)
  
  # (eta - X beta)^\top (D - W - I) (eta - X beta)
  term1 <- as.numeric(crossprod(DWI %*% vector_term, vector_term))
  
  # Tr[(D - W - I) V]
  term2 <- trace_scalars$tr_DWIV
  
  # (term1 + term2) * (-1/2 tau^2)
  constant_term <- (-0.5 / tau2) * (term1 + term2)
  
  # Cleans up computation
  eig_tmp <- eig_vals - 1
  grad_fcn <- function(x) {
    # (1/2) sum_i (kappa_i - 1) / ((kappa_i - 1) gamma + 1)
    eig_term <- 0.5 * sum(eig_tmp / ((eig_tmp * x) + 1.0))
    
    return(eig_term + constant_term)
  }
  
  # Try basic uniroot first
  err_flag <- tryCatch({
    # Find root
    gamma_out <- stats::uniroot(
      grad_fcn,
      lower = 0.0,
      upper = 1.0 - 1e-3,
      # The zero eigenvalue causes issues
      extendInt = "yes"
    )
    
    # Return estimate and gradient value
    return(list(gamma_hat = gamma_out$root, grad_val = gamma_out$f.root))
  }, error = function(cond) {
    warning("UNIROOT FAILED")
    err_flag <- TRUE
  })
  # Then try Newton-Raphson
  if (err_flag) {
    # If null, random
    if (is.null(gamma_current)) {
      gamma_current <- stats::runif(1)
      message("Assigning random gamma ", gamma_current, "\n")
    }
    
    err_flag <- tryCatch({
      gamma_out <- pracma::newtonRaphson(grad_fcn, gamma_current)
      
      # Return estimate and gradient value
      return(list(gamma_hat = gamma_out$root, grad_val = gamma_out$f.root))
    }, error = function(cond) {
      warning("NEWTON-RAPHSON FAILED")
      err_flag <- TRUE
    })
  }
  # Then try fzero
  if (err_flag) {
    # If null, random
    if (is.null(gamma_current)) {
      gamma_current <- stats::runif(1)
      message("Assigning random gamma ", gamma_current, "\n")
    }
    
    err_flag <- tryCatch({
      gamma_out <- pracma::fzero(grad_fcn, gamma_current)
      
      # Return estimate and gradient value
      return(list(gamma_hat = gamma_out$x, grad_val = gamma_out$fval))
    }, error = function(cond) {
      warning("FZERO FAILED")
      err_flag <- TRUE
    })
  }
  if (err_flag) {
    stop("OPTIMIZATION FOR GAMMA FAILED: CANNOT CONTINUE")
  }
}


## Sparse Nearest Neighbor Gaussian Variogram

#' Holding beta constant, fit a variogram model to eta - X beta.
#' Part of the M-Step in the EM algorithm.
#' Instead of a traditional MLE, which has cubic time complexity, we
#' fit a variogram to estimate the kernel parameters, which is quadratic time.
#' ONLY APPLIES TO THE GAUSSIAN PROCESS MODEL.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param eta_hat Estimated mean of eta.
#' @param beta_hat Current estimate of covariate effects.
#' @param X Covariate matrix.
#' @param coords Matrix of (x, y) coordinates.
#' @param cov_type String for covariance model type.
#'  "Exp", "Sph", "Gau", and "Mat" are supported.
#'
#' @returns A vector of fitted parameters.
#'  Nugget variance (partial sill), Spatial variance (partial sill),
#'  Spatial range (scales distance); if Matern, smoothness kappa.
#'
#' @note Requires the gstat library.
#'
#' @import gstat
#' @importFrom gstat variogram
#' @importFrom gstat vgm
#' @importFrom gstat fit.variogram
#' @importFrom stats var
#' @importFrom stats median
M_step_variogram <- function(eta_hat, beta_hat, X, coords, cov_type = "Exp") {
  close_to_zero_const <- max(min(1e-6, 1 / nrow(X)^2), 1e4 * .Machine$double.eps)
  
  # Create a standard base R data frame (bypasses the sp dependency)
  sp_dat <- data.frame(x = coords[, 1],
                       y = coords[, 2],
                       z = as.numeric(eta_hat - X %*% beta_hat))
  
  # Empirical variogram
  # Cutoff is maximum distance: Values are weighted in the fitting
  # Speedup: Vectorized diff(range()) is much faster than apply() or sp::coordinates()
  cutoff <- sqrt(diff(range(sp_dat$x))^2 + diff(range(sp_dat$y))^2)
  
  # Create variogram natively in gstat using the formula's 'locations' argument
  vg_emp <- gstat::variogram(z ~ 1,
                             locations = ~ x + y,
                             data = sp_dat,
                             cutoff = cutoff)
  
  # Define variogram model
  vg_model <- gstat::vgm(
    nugget = min(vg_emp$gamma),
    psill = max(vg_emp$gamma),
    range = stats::median(vg_emp$dist),
    model = cov_type
  )
  
  err_flag <- tryCatch({
    # Fit model
    vg_fit <- gstat::fit.variogram(vg_emp,
                                   vg_model,
                                   fit.kappa = TRUE,
                                   fit.method = 7)
    
    # Nugget variance, Spatial variance, Range, Smoothness (if Matern)
    vg_params <- c(vg_fit[vg_fit$model == "Nug", ]$psill, vg_fit[vg_fit$model == cov_type, ]$psill, vg_fit[vg_fit$model == cov_type, ]$range)
    if ("Mat" == cov_type) {
      vg_params <- c(vg_params, vg_fit[vg_fit$model == cov_type, ]$kappa)
    }
    
    vg_params[is.na(vg_params)] <- close_to_zero_const
    vg_params[is.nan(vg_params)] <- close_to_zero_const
    if ((close_to_zero_const >= vg_params[1]) &&
        (close_to_zero_const >= vg_params[2])) {
      warning(
        "BOTH NUGGET AND SPATIAL VARIANCE ARE ZERO. FAILSAFE: Setting nugget to observed variance."
      )
      vg_params[1] <- max(close_to_zero_const, stats::var(sp_dat$z))
      message("Assigning nugget ", vg_params[1], "\n")
    }
    if (close_to_zero_const >= vg_params[3]) {
      warning("INVALID RANGE; FAILSAFE: Setting range to 95% of p-sill distance.")
      
      # Compute 95% of sill
      sill_95 <- 0.95 * vg_params[2]
      # Find distances corresponding to a value at least 95% of sill
      valid_dists <- which(vg_emp$gamma >= sill_95)
      # If there are valid distances, choose the smallest (first)
      # Otherwise, use the largest distance
      if (1 < length(valid_dists)) {
        range_hat <- vg_emp$dist[valid_dists[1]]
      } else {
        range_hat <- max(vg_emp$dist)
      }
      vg_params[3] <- range_hat
    }
    
    return(vg_params)
  }, error = function(cond) {
    warning("FIT.VARIOGRAM FAILED; WILL DEFAULT TO A FAILSAFE")
    err_flag <- TRUE
  })
  if (err_flag) {
    warning("Initiating FAILSAFE SINCE FIT.VARIOGRAM FAILED: Assuming non-spatial data.")
    vg_params <- c(close_to_zero_const,
                   close_to_zero_const,
                   close_to_zero_const)
    
    # Assigning all variance to non-spatial
    vg_params[1] <- max(close_to_zero_const, stats::var(sp_dat$z))
    
    # Compute 95% of sill
    sill_95 <- 0.95 * vg_params[2]
    # Find distances corresponding to a value at least 95% of sill
    valid_dists <- which(vg_emp$gamma >= sill_95)
    # If there are valid distances, choose the smallest (first)
    # Otherwise, use the largest distance
    if (1 < length(valid_dists)) {
      range_hat <- vg_emp$dist[valid_dists[1]]
    } else {
      range_hat <- max(vg_emp$dist)
    }
    vg_params[3] <- range_hat
    
    # Use a default 3/2 value for kappa
    if ("Mat" == cov_type) {
      vg_params <- c(vg_params, 1.5)
    }
    
    return(vg_params)
  }
}

## Sparse Nearest Neighbor Gaussian Non-Variogram

#' Holding beta constant, use BRISC to find the spatial parameters of eta - X beta.
#' Part of the M-Step in the EM algorithm.
#' ONLY APPLIES TO THE GAUSSIAN PROCESS MODEL.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param eta_hat Estimated mean of eta.
#' @param beta_hat Current estimate of covariate effects.
#' @param X Covariate matrix.
#' @param coords Matrix of (x, y) coordinates.
#' @param cov_type String for covariance model type.
#'  "Exp", "Sph", "Gau", and "Mat" are supported.
#' @param k Number of neighbors.
#'
#' @returns A vector of fitted parameters.
#'  Nugget variance (partial sill), Spatial variance (partial sill),
#'  Spatial range (scales distance); if Matern, smoothness kappa.
#'
#' @note Requires the BRISC library (called indirectly).
#'
#' @importFrom stats var
M_step_BRISC <- function(eta_hat,
                         beta_hat,
                         X,
                         coords,
                         cov_type = "Exp",
                         k = 15) {
  close_to_zero_const <- max(min(1e-6, 1 / nrow(X)^2), 1e4 * .Machine$double.eps)
  
  # Get phi_hat = eta - X beta
  phi_hat <- as.numeric(eta_hat - X %*% beta_hat)
  
  # Call BRISC
  err_flag <- FALSE
  err_flag <- tryCatch({
    b_out <- BRISC_wrapper(
      z_list = list(z = phi_hat),
      # Run on phi
      X_list = NULL,
      # No covariates
      coords_list = list(c = coords),
      # Coordinates
      k = k,
      # Number of neighbors
      cov_type = cov_type,
      # Covariance type
      transform_z = FALSE # Don't log values
    )
    
    if (b_out$Theta[1] + b_out$Theta[2] < close_to_zero_const) {
      warning("Variances are low/unstable; FAILSAFE--Everything is non-spatial")
      b_out$Theta[1] <- max(close_to_zero_const, stats::var(phi_hat))
      b_out$Theta[2] <- max(close_to_zero_const, b_out$Theta[2])
      
      # Use maximum distance as a proxy for the range
      # Speedup: Vectorized diff(range()) instead of apply()
      b_out$Theta[3] <- sqrt(diff(range(coords[, 1]))^2 + diff(range(coords[, 2]))^2)
    }
    return(as.vector(b_out$Theta))
  }, error = function(cond) {
    warning("BRISC FAILED; WILL DEFAULT TO A FAILSAFE")
    err_flag <- TRUE
  })
  if (err_flag) {
    warning("Initiating FAILSAFE SINCE BRISC FAILED: Fitting via a variogram.")
    
    vg_out <- M_step_variogram(eta_hat, beta_hat, X, coords, cov_type)
    return(vg_out)
  }
}

#'  Wrapper to fit sparse NN GP model via BRISC.
#'
#'  This function theoretically works for multiple areas by stacking everything
#'  together and ignoring differences in coordinates.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param z_list List of vectors of observed counts---one vector per area.
#' @param X_list List of design or covariate matrices---one matrix per area.
#'    Same length and ordering as z_list.
#'    Matrices with number of rows equal to length of corresponding vector in
#'    z_list.
#' @param coords_list List of coordinate matrices (x, y)---one matrix per area.
#'    Same length and ordering as z_list.
#'    Matrices with number of rows equal to length of corresponding vector in
#'    z_list.
#' @param k Number of neighbors.
#' @param cov_type String for covariance model type.
#'    "Exp", "Sph", "Gau", and "Mat" are supported.
#' @param transform_z Boolean: log-transform z or not.
#' @param z_offset If transform_z, the counts z are transformed as
#'    log(z + z_offset).
#' @param verbose Whether to print output as BRISC runs.
#'
#' @returns The output list from BRISC, plus a time field for how long the
#'    function ran for.
#'
#' @note Requires the BRISC library.
#'
#' @importFrom BRISC BRISC_estimation
BRISC_wrapper <- function(z_list,
                          X_list,
                          coords_list,
                          k = 15,
                          cov_type = "Exp",
                          transform_z = TRUE,
                          z_offset = 0.5,
                          verbose = FALSE) {
  # Start clock
  t0_glm <- Sys.time()
  
  # Stack into a single vector/matrix efficiently
  # Speedup: unlist() and do.call(rbind, ...) bypass Reduce() memory bloat
  z_vec <- unlist(z_list, use.names = FALSE)
  X_mat <- if (!is.null(X_list))
    do.call(rbind, X_list)
  else
    NULL
  coords <- do.call(rbind, coords_list)
  
  # Transform z
  if (transform_z) {
    z_vec <- log(z_vec + z_offset)
  }
  
  if ("Exp" == cov_type) {
    brisc_cov <- "exponential"
  } else if ("Sph" == cov_type) {
    brisc_cov <- "spherical"
  } else if ("Mat" == cov_type) {
    brisc_cov <- "matern"
  } else if ("Gau" == cov_type) {
    brisc_cov <- "gaussian"
  }
  
  # Fit LM
  b_out <- BRISC::BRISC_estimation(
    coords,
    z_vec,
    X_mat,
    n.neighbors = k,
    cov.model = brisc_cov,
    verbose = verbose
  )
  
  # Stop clock
  t1_glm <- Sys.time()
  b_out[["time"]] <- difftime(t1_glm, t0_glm)
  
  return(b_out)
}
