## M-Step (Maximization Step) functions.
# Dependencies in file: Matrix, pracma, stats, sp, gstat.
# Dependencies: Functions from comparisons.R.

#' Maximize the expected likelihood in tau^2, holding other variables constant.
#' Part of the M-Step in the EM algorithm.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param Vhat Estimated covariance matrix of eta.
#' @param eta_hat Estimated mean of eta.
#' @param Q Unscaled precision matrix.
#' @param beta_hat Current estimate of covariate effects.
#' @param X Covariate matrix.
#'
#' @returns Maximum likelihood estimate of scaling parameter tau^2.
M_step_tau2 <- function(Vhat, eta_hat, Q, beta_hat, X) {
  # eta - X beta
  vector_term <- eta_hat - (X %*% beta_hat)

  # (eta - X beta)^\top Q (eta - X beta)
  # term1 <- as.numeric((t(vector_term) %*% Q) %*% vector_term)
  # Speed
  term1 <- as.numeric(crossprod(Q %*% vector_term, vector_term))

  # Trace[Q V]
  # term2 <- sum(diag(Q %*% Vhat))
  # Trace tricks: Tr(A B) = <vec(A), vec(B^T)>
  # For sparse Q this is actually faster than crossprod?
  term2 <- sum(Q * t(Vhat))
  # term2 <- crossprod(as.vector(Q), as.vector(t(Vhat)))

  # tau^2 = (1/n) (term1 + term2)
  tau2_hat <- (term1 + term2) / dim(Q)[1]

  # Make sure tau^2 doesn't go off the rails
  close_to_zero_const <- 1e2 * .Machine$double.eps
  if (close_to_zero_const >= tau2_hat) {
    print("Invalid tau^2, setting to a small positive number")
    tau2_hat = max(tau2_hat, close_to_zero_const)
  }

  return(tau2_hat)
}

#' Maximize the expected likelihood in beta, holding other variables constant.
#' Part of the M-Step in the EM algorithm.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to ALL areas.
#'
#' @param eta_list List of the estimated means of eta.
#'  One vector per area.
#' @param Q_list List of the unscaled precision matrices.
#'  One matrix per area, same ordering and length as eta_list.
#'  Depends on gamma, so that dependency is implicit.
#'  I.e., Q_list should be updated in each iteration after gamma is updated.
#' @param tau2_list List or vector of precision matrix scaling values.
#'  Same length and ordering as eta_list, with one number per area.
#' @param X_list List of covariate matrices.
#'  One matrix per area, same ordering and length as eta_list.
#'
#' @returns Maximum likelihood estimate of covariates beta.
#'
#' @note Requires the Matrix library.
#' @note Requires the pracma library.
#'
#' @import Matrix
#' @importFrom Matrix solve
#' @importFrom pracma pinv
M_step_beta <- function(eta_list, Q_list, tau2_list, X_list) {
  # Dimension of beta
  n_dim <- ncol(X_list[[1]])

  # Create empty matrix/vector
  zeta_vec <- matrix(0, nrow = n_dim, ncol = 1)
  B <- matrix(0, nrow = n_dim, ncol = n_dim)
  for (idx in 1:length(Q_list)) {
    # X^\top Q
    # term1 <- t(X_list[[idx]]) %*% Q_list[[idx]]
    term1 <- crossprod(X_list[[idx]], Q_list[[idx]])

    # Update vector term
    # X^\top (Q / tau^2) eta = (X^\top Q) eta / tau^2
    zeta_vec <-
      zeta_vec + (term1 %*% eta_list[[idx]]) / tau2_list[[idx]]

    # Update matrix term
    # X^\top (Q / tau^2) X = (X^\top Q) X / tau^2
    B <- B + (term1 %*% X_list[[idx]]) / tau2_list[[idx]]
  }

  # Try basic inversion first
  err_flag <- tryCatch({
    beta_hat <- Matrix::solve(B, zeta_vec)
    return(beta_hat)
  }, error = function(cond) {
    print("INVERSION OF B FAILED")
    err_flag <- TRUE
  })
  # Then try pseudoinverse
  if (err_flag) {
    err_flag <- tryCatch({
      beta_hat <- pracma::pinv(as.matrix(B)) %*% zeta_vec
      return(beta_hat)
    }, error = function(cond) {
      print("PSEUDOINVERSE OF B FAILED")
      err_flag <- TRUE
    })
  }
  if (err_flag) {
    stop("PSEUDOINVERSE OF B FAILED: CANNOT CONTINUE")
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
#' @param Vhat Estimated covariance matrix of eta.
#' @param eta_hat Estimated mean of eta.
#' @param tau2 Precision matrix scaling parameter tau^2.
#' @param beta_hat Current estimate of covariate effects.
#' @param X Covariate matrix.
#' @param W Neighbor/adjacency matrix (symmetric, binary).
#' @param D Diagonal degree matrix (row sums of W).
#'  Not needed, but passed in for consistency with API.
#' @param eig_vals Eigenvalues of Z = D^{-1} W.
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
M_step_gamma_CAR <- function(Vhat,
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
  # term1 <- as.numeric((t(vector_term) %*% W) %*% vector_term)
  # Speed
  term1 <- as.numeric(crossprod(W %*% vector_term, vector_term))

  # Tr[W V]
  # term2 <- sum(diag(W %*% Vhat))
  # Trace tricks
  # term2 <- crossprod(as.vector(W), as.vector(t(Vhat)))
  # For sparse W this is actually faster
  term2 <- sum(W * t(Vhat))

  # (term1 + term2) / (2 tau^2)
  constant_term <- (term1 + term2) * (0.5 / tau2)

  # Gradient function
  # Inner function that evaluates the gradient
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
    print("UNIROOT FAILED")
    err_flag <- TRUE
  })
  # Then try Newton-Raphson
  if (err_flag) {
    # If null, random
    if (is.null(gamma_current)) {
      gamma_current = stats::runif(1)
      print(gamma_current)
    }

    err_flag <- tryCatch({
      gamma_out <- pracma::newtonRaphson(grad_fcn, gamma_current)

      # print(gamma_out) # Since it failed, let's diagnose it

      # Return estimate and gradient value
      return(list(gamma_hat = gamma_out$root, grad_val = gamma_out$f.root))
    }, error = function(cond) {
      print("NEWTON-RAPHSON FAILED")
      err_flag <- TRUE
    })
  }
  # Then try fzero
  if (err_flag) {
    # If null, random
    if (is.null(gamma_current)) {
      gamma_current = stats::runif(1)
      print(gamma_current)
    }

    err_flag <- tryCatch({
      gamma_out <- pracma::fzero(grad_fcn, gamma_current)

      # print(gamma_out) # Since it failed, let's diagnose it

      # Return estimate and gradient value
      return(list(gamma_hat = gamma_out$x, grad_val = gamma_out$fval))
    }, error = function(cond) {
      print("FZERO FAILED")
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
#' @param Vhat Estimated covariance matrix of eta.
#' @param eta_hat Estimated mean of eta.
#' @param tau2 Precision matrix scaling parameter tau^2.
#' @param beta_hat Current estimate of covariate effects.
#' @param X Covariate matrix.
#' @param W Neighbor/adjacency matrix (symmetric, binary).
#' @param D Diagonal degree matrix (row sums of W).
#' @param eig_vals Eigenvalues of Z = D^{-1} W.
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
M_step_gamma_SAR <- function(Vhat,
                             eta_hat,
                             tau2,
                             beta_hat,
                             X,
                             W,
                             D,
                             eig_vals,
                             gamma_current = NULL) {
  # D^{-1}
  D_inv <-
    Matrix::Diagonal(dim(W)[1], 1 / Matrix::diag(D))
  # Z = D^{-1} W
  Z <- D_inv %*% W

  # eta - X beta
  vector_term <- eta_hat - (X %*% beta_hat)

  # (eta - X beta)^\top (2 W) (eta - X beta)
  # zeta1 <- 2.0 * as.numeric((t(vector_term) %*% W) %*% vector_term)
  # Speed
  zeta1 <- 2.0 * as.numeric(crossprod(W %*% vector_term, vector_term))

  # (eta - X beta)^\top W Z (eta - X beta)
  # zeta2 <- as.numeric((t(vector_term) %*% W) %*% (Z %*% vector_term))
  zeta2 <- as.numeric(crossprod(W %*% vector_term, Z %*% vector_term))

  # Tr[W V]
  # zeta3 <- sum(diag(W %*% Vhat))
  # Trace tricks
  # zeta3 <- crossprod(as.vector(W), as.vector(t(Vhat)))
  # This is faster for sparse W
  zeta3 <- sum(W * t(Vhat))

  # Tr[W Z V]
  # zeta4 <- sum(diag((W %*% Z) %*% Vhat))
  # Trace tricks: Tr(A B) = <vec(A), vec(B^T)>
  # Also note that Z V becomes sparse, so this grouping is faster
  # This is faster for sparse W and Z
  zeta4 <- sum(W * t(Z %*% Vhat))
  # An even faster method
  # zeta4 <- crossprod(as.vector(W), as.vector(t(Z %*% Vhat)))

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
    print("UNIROOT FAILED")
    err_flag <- TRUE
  })
  # Then try Newton-Raphson
  if (err_flag) {
    # If null, random
    if (is.null(gamma_current)) {
      gamma_current = stats::runif(1)
      print(gamma_current)
    }

    err_flag <- tryCatch({
      gamma_out <- pracma::newtonRaphson(grad_fcn, gamma_current)

      # print(gamma_out) # Since it failed, let's diagnose it

      # Return estimate and gradient value
      return(list(gamma_hat = gamma_out$root, grad_val = gamma_out$f.root))
    }, error = function(cond) {
      print("NEWTON-RAPHSON FAILED")
      err_flag <- TRUE
    })
  }
  # Then try fzero
  if (err_flag) {
    # If null, random
    if (is.null(gamma_current)) {
      gamma_current = stats::runif(1)
      print(gamma_current)
    }

    err_flag <- tryCatch({
      gamma_out <- pracma::fzero(grad_fcn, gamma_current)

      # print(gamma_out) # Since it failed, let's diagnose it

      # Return estimate and gradient value
      return(list(gamma_hat = gamma_out$x, grad_val = gamma_out$fval))
    }, error = function(cond) {
      print("FZERO FAILED")
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
#' @param Vhat Estimated covariance matrix of eta.
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
M_step_gamma_Leroux <- function(Vhat,
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
  # term1 <- as.numeric((t(vector_term) %*% DWI) %*% vector_term)
  # Speed
  term1 <- as.numeric(crossprod(DWI %*% vector_term, vector_term))

  # Tr[(D - W - I) V]
  # term2 <- sum(diag(DWI %*% Vhat))
  # Trace tricks
  # term2 <- as.numeric(crossprod(as.vector(DWI), as.vector(t(Vhat))))
  # This is faster for sparse D and W
  term2 <- sum(DWI * t(Vhat))

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
    print("UNIROOT FAILED")
    err_flag <- TRUE
  })
  # Then try Newton-Raphson
  if (err_flag) {
    # If null, random
    if (is.null(gamma_current)) {
      gamma_current = stats::runif(1)
      print(gamma_current)
    }

    err_flag <- tryCatch({
      gamma_out <- pracma::newtonRaphson(grad_fcn, gamma_current)

      # print(gamma_out) # Since it failed, let's diagnose it

      # Return estimate and gradient value
      return(list(gamma_hat = gamma_out$root, grad_val = gamma_out$f.root))
    }, error = function(cond) {
      print("NEWTON-RAPHSON FAILED")
      err_flag <- TRUE
    })
  }
  # Then try fzero
  if (err_flag) {
    # If null, random
    if (is.null(gamma_current)) {
      gamma_current = stats::runif(1)
      print(gamma_current)
    }

    err_flag <- tryCatch({
      gamma_out <- pracma::fzero(grad_fcn, gamma_current)

      # print(gamma_out) # Since it failed, let's diagnose it

      # Return estimate and gradient value
      return(list(gamma_hat = gamma_out$x, grad_val = gamma_out$fval))
    }, error = function(cond) {
      print("FZERO FAILED")
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
#' @note Requires the sp library.
#' @note Requires the gstat library.
#'
#' @import sp
#' @importFrom sp coordinates
#' @import gstat
#' @importFrom gstat variogram
#' @importFrom gstat vgm
#' @importFrom gstat fit.variogram
#' @importFrom stats var
#' @importFrom stats median
M_step_variogram <- function(eta_hat, beta_hat, X, coords, cov_type = "Exp") {
  # close_to_zero_const <- 1e4 * .Machine$double.eps
  close_to_zero_const <- max(min(1e-6, 1 / nrow(X)^2), 1e4 * .Machine$double.eps)

  # Create data frame for spatial packages
  sp_dat <- cbind(coords, as.numeric(eta_hat - X %*% beta_hat))
  colnames(sp_dat) <- c("x", "y", "z")
  sp_dat <- as.data.frame(sp_dat)

  # Identify coordinates for spatial packages
  sp::coordinates(sp_dat) <- ~ x + y

  # Empirical variogram
  # Cutoff is maximum distance: Values are weighted in the fitting
  cutoff <- sqrt(diff(range(sp::coordinates(sp_dat)[, 1]))^2
                 + diff(range(sp::coordinates(sp_dat)[, 2]))^2)
  # Create variogram
  vg_emp <- gstat::variogram(z ~ 1, data = sp_dat, cutoff = cutoff)

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
      print("BOTH NUGGET AND SPATIAL VARIANCE ARE ZERO.")
      print("FAILSAFE: Setting nugget to observed variance.")
      vg_params[1] <- max(close_to_zero_const, stats::var(sp_dat$z))
      print(vg_params[1])
    }
    if (close_to_zero_const >= vg_params[3]) {
      print("WARNING: INVALID RANGE")
      print("FAILSAFE: Setting range to 95% of p-sill distance.")

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
    print("FIT.VARIOGRAM FAILED")
    print("WILL DEFAULT TO A FAILSAFE")
    err_flag <- TRUE
  })
  if (err_flag) {
    print("Initiating FAILSAFE SINCE FIT.VARIOGRAM FAILED: Assuming non-spatial data.")
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
  # close_to_zero_const <- 1e4 * .Machine$double.eps
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
      print("Variances are low/unstable; FAILSAFE--Everything is non-spatial")
      b_out$Theta[1] <- max(close_to_zero_const, stats::var(phi_hat))
      b_out$Theta[2] <- max(close_to_zero_const, b_out$Theta[2])
      # Use maximum distance as a proxy for the range
      b_out$Theta[3] <- sqrt(sum((
        apply(coords, 2, max) - apply(coords, 2, min)
      )^2))
    }
    return(as.vector(b_out$Theta))
  }, error = function(cond) {
    print("BRISC FAILED")
    print("WILL DEFAULT TO A FAILSAFE")
    err_flag <- TRUE
  })
  if (err_flag) {
    print("Initiating FAILSAFE SINCE BRISC FAILED: Fitting via a variogram.")

    # Return zeros for the variances, and default values for the rest
    # close_to_zero_const <- 2.0 * .Machine$double.eps
    # return( c(max(close_to_zero_const, stats::var(phi_hat)), close_to_zero_const, 1, 1.5) )

    vg_out <- M_step_variogram(eta_hat, beta_hat, X, coords, cov_type)
    return(vg_out)
  }
}
