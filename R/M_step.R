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


#' Maximize the expected likelihood in beta, holding other variables constant.
#' Part of the M-Step in the EM algorithm.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies across all areas.
#'
#' @param eta_list List of estimated means of eta across all areas.
#' @param Q_list List of unscaled precision matrices across all areas.
#' @param tau2_list Vector of precision matrix scaling parameters tau^2 across all areas.
#' @param X_list List of covariate matrices across all areas.
#' @param model_type Character string specifying "CAR", "SAR", "Leroux", or "spNNGP".
#' @param gamma_list Vector of current spatial correlation parameters gamma across all areas.
#' @param XtDX_list Optional list of precomputed X^T D X matrices (CAR/SAR).
#' @param XtWX_list Optional list of precomputed X^T W X matrices (CAR/SAR).
#' @param XtWZX_list Optional list of precomputed X^T W Z X matrices (SAR).
#' @param XtX_list Optional list of precomputed X^T X matrices (Leroux).
#' @param XtDWIX_list Optional list of precomputed X^T (D - W - I) X matrices (Leroux).
#' @param XtQX_list Optional list of precomputed X^T Q X matrices (spNNGP inner CM steps).
#' @param XtQeta_list Optional list of precomputed X^T Q eta vectors (spNNGP inner CM steps).
#'
#' @returns Maximum likelihood estimate of covariate effects beta.
#'
#' @note Requires the Matrix library.
#' @note Requires the pracma library.
#'
#' @import Matrix
#' @importFrom Matrix crossprod
#' @importFrom pracma pinv
M_step_beta <- function(eta_list,
                        Q_list,
                        tau2_list,
                        X_list,
                        model_type = "Leroux",
                        gamma_list = NULL,
                        XtDX_list = NULL,
                        XtWX_list = NULL,
                        XtWZX_list = NULL,
                        XtX_list = NULL,
                        XtDWIX_list = NULL,
                        XtQX_list = NULL,
                        XtQeta_list = NULL) {
  # Get number of covariates
  p <- ncol(X_list[[1]])
  # Will store sums of matrices/vectors across areas
  B_mat <- matrix(0, nrow = p, ncol = p)
  zeta_val <- rep(0, p)
  
  # Iterate over areas
  for (idx in 1:length(eta_list)) {
    X <- X_list[[idx]]
    eta <- eta_list[[idx]]
    tau2 <- tau2_list[idx]
    Q <- Q_list[[idx]]
    gamma <- gamma_list[idx]
    
    # Fast path for spNNGP when X^T Q eta and X^T Q X were computed directly in C++
    if (model_type == "spNNGP" && !is.null(XtQX_list) && !is.null(XtQeta_list)) {
      zeta_val <- zeta_val + (XtQeta_list[[idx]] / tau2)
      B_mat <- B_mat + (XtQX_list[[idx]] / tau2)
      next
    }
    
    # 1. Compute RHS: X^T Q eta / tau^2
    # Optimized to avoid forming dense intermediate products; Q %*% eta is O(N)
    zeta_val <- zeta_val + (Matrix::crossprod(X, Q %*% eta)) / tau2
    
    # 2. Compute LHS: X^T Q X / tau^2
    # Utilizing precomputed static matrices if provided for the specific lattice model
    if (model_type == "CAR" &&
        !is.null(XtDX_list) && !is.null(XtWX_list)) {
      # Q = D - gamma * W
      XtQX <- XtDX_list[[idx]] - gamma * XtWX_list[[idx]]
      B_mat <- B_mat + (XtQX / tau2)
      
    } else if (model_type == "SAR" &&
               !is.null(XtDX_list) &&
               !is.null(XtWX_list) && !is.null(XtWZX_list)) {
      # Q = D - 2*gamma*W + gamma^2 * W * D^{-1} * W
      XtQX <- XtDX_list[[idx]] - (2.0 * gamma) * XtWX_list[[idx]] + (gamma^2) * XtWZX_list[[idx]]
      B_mat <- B_mat + (XtQX / tau2)
      
    } else if (model_type == "Leroux" &&
               !is.null(XtX_list) && !is.null(XtDWIX_list)) {
      # Q = I + gamma * (D - W - I)
      XtQX <- XtX_list[[idx]] + gamma * XtDWIX_list[[idx]]
      B_mat <- B_mat + (XtQX / tau2)
      
    } else {
      # Fallback for spNNGP (when Q is explicitly formed) or if lists are missing
      B_mat <- B_mat + (Matrix::crossprod(X, Q %*% X)) / tau2
    }
  }
  
  # Ensure B_mat is a standard base matrix for fast inversion (p x p is small)
  B_mat <- as.matrix(B_mat)
  zeta_val <- as.numeric(zeta_val)
  
  # Try basic inversion first
  err_flag <- FALSE
  err_flag <- tryCatch({
    beta_out <- solve(B_mat, zeta_val)
    return(as.numeric(beta_out))
  }, error = function(cond) {
    warning("INVERSION OF B FAILED")
    err_flag <- TRUE
  })
  
  # Try pseudoinversion next
  if (err_flag) {
    err_flag <- tryCatch({
      B_inv <- pracma::pinv(B_mat)
      return(as.numeric(B_inv %*% zeta_val))
    }, error = function(cond) {
      warning("PSEUDOINVERSION OF B FAILED")
      err_flag <- TRUE
    })
  }
  if (err_flag) {
    stop("ALL INVERSIONS OF B FAILED: CANNOT CONTINUE")
  }
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
#' @param Z_precomp Optional precomputed Z = D^\{-1\} W matrix.
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
                             gamma_current = NULL,
                             Z_precomp = NULL) {
  if (is.null(Z_precomp)) {
    # D^{-1}
    D_inv <- Matrix::Diagonal(dim(W)[1], 1 / Matrix::diag(D))
    # Z = D^{-1} W
    Z <- D_inv %*% W
  } else {
    Z <- Z_precomp
  }
  
  # eta - X beta
  vector_term <- eta_hat - (X %*% beta_hat)
  
  # (eta - X beta)^\top (2 W) (eta - X beta)
  Wg <- W %*% vector_term
  zeta1 <- 2.0 * as.numeric(crossprod(Wg, vector_term))
  
  # (eta - X beta)^\top W Z (eta - X beta)
  zeta2 <- as.numeric(crossprod(Wg, Z %*% vector_term))
  
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
#' @param DWI_precomp Optional precomputed D - W - I matrix.
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
                                gamma_current = NULL,
                                DWI_precomp = NULL) {
  if (is.null(DWI_precomp)) {
    # Identity matrix
    id_mat <- Matrix::Diagonal(dim(W)[1], 1)
    # D - W - I
    DWI <- D - W - id_mat
  } else {
    DWI <- DWI_precomp
  }
  
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
#' @param prev_params Optional numeric vector of previous covariance parameters
#'  in gstat convention: c(nugget, sill, range, \[smoothness\]). Used to warm-start
#'  BRISC optimization across EM/CM iterations.
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
                         k = 15,
                         prev_params = NULL) {
  close_to_zero_const <- max(min(1e-6, 1 / nrow(X)^2), 1e4 * .Machine$double.eps)
  
  # Get phi_hat = eta - X beta
  phi_hat <- as.numeric(eta_hat - X %*% beta_hat)
  
  # Call BRISC (try warm-start first if prev_params supplied; fallback to cold-start)
  err_flag <- FALSE
  err_flag <- tryCatch({
    b_out <- tryCatch({
      BRISC_wrapper(
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
        transform_z = FALSE,
        # Don't log values
        prev_params = prev_params
      )
    }, error = function(cond_warm) {
      # If warm-start fails, retry with BRISC default cold-start before variogram fallback
      if (!is.null(prev_params)) {
        BRISC_wrapper(
          z_list = list(z = phi_hat),
          X_list = NULL,
          coords_list = list(c = coords),
          k = k,
          cov_type = cov_type,
          transform_z = FALSE,
          prev_params = NULL
        )
      } else {
        stop(cond_warm)
      }
    })
    
    # Initialize a new vector for gstat convention: c(Nugget, Sill, Range, [Smoothness])
    out_params <- numeric(length(b_out$Theta))
    
    # Check if total variance is near zero using BRISC's native indices (1=Sill, 2=Nugget)
    if (b_out$Theta[1] + b_out$Theta[2] < close_to_zero_const) {
      warning("Variances are low/unstable; FAILSAFE--Everything is non-spatial")
      
      # gstat order: Index 1 = Nugget, Index 2 = Sill, Index 3 = Range
      out_params[1] <- max(close_to_zero_const, b_out$Theta[2])
      out_params[2] <- max(close_to_zero_const, stats::var(phi_hat))
      
      # Use maximum distance as a proxy for the range
      # Speedup: Vectorized diff(range()) instead of apply()
      out_params[3] <- sqrt(diff(range(coords[, 1]))^2 + diff(range(coords[, 2]))^2)
      
    } else {
      # Map BRISC output to gstat convention
      out_params[1] <- b_out$Theta[2]       # Nugget (BRISC tau.sq)
      out_params[2] <- b_out$Theta[1]       # Sill (BRISC sigma.sq)
      out_params[3] <- 1.0 / b_out$Theta[3] # Range (1 / BRISC phi)
    }
    
    # If a 4th parameter exists (Smoothness / nu for Matern), pass it through
    if (length(b_out$Theta) == 4) {
      out_params[4] <- b_out$Theta[4]
    }
    
    return(out_params)
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

#' Wrapper to fit sparse NN GP model via BRISC.
#'
#' This function theoretically works for multiple areas by stacking everything
#' together and ignoring differences in coordinates.
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
#' @param prev_params Optional numeric vector of previous covariance parameters
#'    in gstat convention: c(nugget, sill, range, \[smoothness\]). When supplied,
#'    maps to BRISC's (sigma.sq, tau.sq, phi, nu) starting values.
#'
#' @returns The output list from BRISC, plus a time field for how long the
#'    function ran for.
#'
#' @note Requires the BRISC library.
#' @note BRISC uses a different convention than gstat for ordering and defining
#'  its spatial covariance parameters: gstat uses (nugget, sill, range, order),
#'  but BRISC uses (sill, nugget, 1 / range, order). This function does NOT change
#'  that, but we (internally) make adjustments where we call it.
#'
#' @importFrom BRISC BRISC_estimation
BRISC_wrapper <- function(z_list,
                          X_list,
                          coords_list,
                          k = 15,
                          cov_type = "Exp",
                          transform_z = TRUE,
                          z_offset = 0.5,
                          verbose = FALSE,
                          prev_params = NULL) {
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
  
  # Default BRISC starting values
  sigma_sq_init <- 1.0
  tau_sq_init   <- 0.1
  phi_init      <- 1.0
  nu_init       <- 1.5
  
  # Override with warm-start values if valid previous gstat parameters exist:
  # prev_params[1] = Nugget (tau.sq), prev_params[2] = Sill (sigma.sq),
  # prev_params[3] = Range (1 / phi), prev_params[4] = Smoothness (nu)
  if (!is.null(prev_params) && length(prev_params) >= 3 &&
      all(is.finite(prev_params[1:3])) &&
      all(prev_params[1:3] > 0)) {
    tau_sq_init   <- max(1e-6, prev_params[1])
    sigma_sq_init <- max(1e-6, prev_params[2])
    phi_init      <- max(1e-6, 1.0 / prev_params[3])
    if ("Mat" == cov_type && length(prev_params) >= 4 &&
        is.finite(prev_params[4]) && prev_params[4] > 0) {
      nu_init <- min(max(1e-3, prev_params[4]), 10.0)
    }
  }
  
  # Fit LM via BRISC (positional args 1-3 match coords, y, x)
  b_out <- BRISC::BRISC_estimation(
    coords,
    z_vec,
    X_mat,
    sigma.sq = sigma_sq_init,
    tau.sq = tau_sq_init,
    phi = phi_init,
    nu = nu_init,
    n.neighbors = k,
    cov.model = brisc_cov,
    verbose = verbose
  )
  
  # Stop clock
  t1_glm <- Sys.time()
  b_out[["time"]] <- difftime(t1_glm, t0_glm)
  
  return(b_out)
}
