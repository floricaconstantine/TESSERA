## Model-specific functions, e.g., precision matrices.
# Dependencies in file: Matrix.

## Precision matrix functions.

#' Compute the unscaled precision matrix in a CAR model.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param W Neighbor/adjacency matrix (symmetric, binary).
#' @param D Degree matrix (diagonal, values are row-sums of W).
#' @param gamma_val Correlation parameter.
#'
#' @returns Unscaled precision matrix.
#' @export
Q_matrix_CAR <- function(W, D, gamma_val) {
  return(D - gamma_val * W)
}

#' Compute the unscaled precision matrix in a SAR model.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param W Neighbor/adjacency matrix (symmetric, binary).
#' @param D Degree matrix (diagonal, values are row-sums of W).
#' @param gamma_val Correlation parameter.
#'
#' @returns Unscaled precision matrix.
#'
#' @note Requires the Matrix library.
#'
#' @import Matrix
#' @importFrom Matrix Diagonal
#' @importFrom Matrix diag
#'
#' @export
Q_matrix_SAR <- function(W, D, gamma_val) {
  # D^{-1}
  D_inv <-
    Matrix::Diagonal(dim(W)[1], 1 / Matrix::diag(D))
  # Identity matrix
  id_mat <- Matrix::Diagonal(dim(W)[1], 1)
  # I - gamma Z = I - gamma D^{-1} W
  I_gamma_Z <- id_mat - ((gamma_val * D_inv) %*% W)

  # (I - gamma Z)^\top D (I - gamma Z)
  return((t(I_gamma_Z) %*% D) %*% I_gamma_Z)
}

#' Compute the unscaled precision matrix in a Leroux model.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param W Neighbor/adjacency matrix (symmetric, binary).
#' @param D Degree matrix (diagonal, values are row-sums of W).
#' @param gamma_val Correlation parameter.
#'
#' @returns Unscaled precision matrix.
#'
#' @note Requires the Matrix library.
#'
#' @import Matrix
#' @importFrom Matrix Diagonal
#'
#' @export
Q_matrix_Leroux <- function(W, D, gamma_val) {
  # Identity matrix
  id_mat <- Matrix::Diagonal(dim(W)[1], 1)

  # gamma (D - W) + (1 - gamma) I
  return(gamma_val * (D - W) + (1 - gamma_val) * id_mat)
}

## Negative Hessians (for Wald tests, etc.)

#' Compute the Negative Hessian for beta.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to ALL areas.
#'
#' @param Q_list List of the unscaled precision matrices.
#'  One matrix per area, same ordering and length as eta_list.
#'  Depends on gamma, so that dependency is implicit.
#'  I.e., Q_list should be updated in each iteration after gamma is updated.
#' @param tau2_list List or vector of precision matrix scaling values.
#'  Same length and ordering as eta_list, with one number per area.
#' @param X_list List of covariate matrices.
#'  One matrix per area, same ordering and length as eta_list.
#'
#' @returns Negative Hessian matrix of beta.
#'  Also called the observed information matrix.
#'
#' @note Requires the Matrix library.
#'
#' @import Matrix
neg_hessian_beta <- function(Q_list, tau2_list, X_list) {
  # Dimension of beta
  n_dim <- ncol(X_list[[1]])

  # Create empty matrix
  B <- matrix(0, nrow = n_dim, ncol = n_dim)
  for (idx in 1:length(Q_list)) {
    # Update matrix term
    # X^\top (Q / tau^2) X = (X^\top Q) X / tau^2
    B <- B + crossprod(X_list[[idx]], Q_list[[idx]] %*% X_list[[idx]]) / tau2_list[[idx]]
  }
  return(B)
}

#' Compute the Negative Hessian for tau^2.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param Vhat Estimated covariance matrix of eta.
#' @param eta_hat Estimated mean of eta.
#' @param Q Unscaled precision matrix.
#' @param tau2_hat Parameter of interest: Current estimate of scale.
#' @param beta_hat Current estimate of covariate effects.
#' @param X Covariate matrix.
#'
#' @returns Negative Hessian of tau^2.
#'  Also called the observed information.
neg_hessian_tau2 <- function(Vhat, eta_hat, Q, tau2_hat, beta_hat, X) {
  # eta - X beta
  vector_term <- eta_hat - (X %*% beta_hat)

  # (eta - X beta)^\top Q (eta - X beta)
  # term1 <- as.numeric((t(vector_term) %*% Q) %*% vector_term)
  # Speed
  term1 <- as.numeric(crossprod(Q %*% vector_term, vector_term))

  # Trace[Q V]
  # term2 <- sum(diag(Q %*% Vhat))
  # Trace tricks: Tr(A B) = <vec(A), vec(B^T)>
  # Faster for sparse Q
  term2 <- sum(Q * t(Vhat))
  # term2 <- crossprod(as.vector(Q), as.vector(t(Vhat)))

  return((term1 + term2) / (tau2_hat^3) - (0.5 * nrow(X)) / (tau2_hat^2))
}

#' Compute the Negative Hessian of gamma.
#' ONLY APPLIES TO THE CAR MODEL.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param gamma_hat Parameter of interest: Estimate of correlation parameter.
#' @param eig_vals Eigenvalues of Z = D^{-1} W.
#'
#' @returns Negative Hessian of gamma.
#'  Also called the observed information.
neg_hessian_gamma_CAR <- function(gamma_hat, eig_vals) {
  return(0.5 * sum((eig_vals / (
    1.0 - (gamma_hat * eig_vals)
  ))^2))
}

#' Compute the Negative Hessian of gamma.
#' ONLY APPLIES TO THE SAR MODEL.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param Vhat Estimated covariance matrix of eta.
#' @param eta_hat Estimated mean of eta.
#' @param gamma_hat Parameter of interest: Estimate of correlation parameter.
#' @param tau2 Parameter of interest: Current estimate of scale.
#' @param beta_hat Current estimate of covariate effects.
#' @param X Covariate matrix.
#' @param W Neighbor/adjacency matrix (symmetric, binary).
#' @param D Degree matrix (diagonal, values are row-sums of W).
#' @param eig_vals Eigenvalues of Z = D^{-1} W.
#'
#' @returns Negative Hessian of gamma.
#'  Also called the observed information.
#'
#' @note Requires the Matrix library.
#'
#' @import Matrix
#' @importFrom Matrix Diagonal
#' @importFrom Matrix diag
neg_hessian_gamma_SAR <- function(Vhat,
                                  eta_hat,
                                  gamma_hat,
                                  tau2,
                                  beta_hat,
                                  X,
                                  W,
                                  D,
                                  eig_vals) {
  eig_term <- sum((eig_vals / (1.0 - (
    gamma_hat * eig_vals
  )))^2)

  # D^{-1}
  D_inv <-
    Matrix::Diagonal(dim(W)[1], 1 / Matrix::diag(D))
  # Z = D^{-1} W
  Z <- D_inv %*% W

  # eta - X beta
  vector_term <- eta_hat - (X %*% beta_hat)

  # (eta - X beta)^\top W Z (eta - X beta)
  # zeta2 <- as.numeric((t(vector_term) %*% W) %*% (Z %*% vector_term))
  zeta2 <- as.numeric(crossprod(W %*% vector_term, Z %*% vector_term))

  # Tr[W Z V]
  # zeta4 <- sum(diag((W %*% Z) %*% Vhat))
  # Trace tricks: Tr(A B) = <vec(A), vec(B^T)>
  # Also note that Z V becomes sparse, so this grouping is faster
  # This is faster for sparse W and Z
  zeta4 <- sum(W * t(Z %*% Vhat))
  # An even faster method
  # zeta4 <- crossprod(as.vector(W), as.vector(t(Z %*% Vhat)))

  return(eig_term + (zeta2 + zeta4) / tau2)
}

#' Compute the Negative Hessian of gamma.
#' ONLY APPLIES TO THE Leroux MODEL.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param gamma_hat Parameter of interest: Estimate of correlation parameter.
#' @param eig_vals Eigenvalues of D - W.
#'
#' @returns Negative Hessian of gamma.
#'  Also called the observed information.
neg_hessian_gamma_Leroux <- function(gamma_hat, eig_vals) {
  eig_tmp <- eig_vals - 1
  return(sum((eig_tmp / (
    eig_tmp * gamma_hat + 1
  ))^2))
}

## Expected log likelihood

#' Compute the expected log likelihood.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#' @note Multivariate normal; drops the 2 pi term.
#'
#' @param Vhat Estimated covariance matrix of eta.
#' @param eta_hat Estimated mean of eta.
#' @param Q Unscaled precision matrix.
#'  The scaled precision matrix in spNNGP.
#' @param gamma_hat Parameter of interest: Estimate of correlation parameter.
#' @param tau2 Parameter of interest: Current estimate of scale.
#' @param beta_hat Current estimate of covariate effects.
#' @param X Covariate matrix.
#' @param W Neighbor/adjacency matrix (symmetric, binary).
#' @param D Degree matrix (diagonal, values are row-sums of W).
#' @param eig_vals Eigenvalues of Z = D^{-1} W (CAR/SAR), D - W (Leroux), or Q (spNNGP).
#' @param model_type "CAR", "SAR", or "Leroux", or "spNNGP".
#'    If spNNGP, set tau^2 to 1 and the value of gamma is irrelevant.
#'  Model for the random effects.
#'
#' @returns Expected log likelihood.
expected_loglike <- function(Vhat,
                             eta_hat,
                             Q,
                             gamma_hat,
                             tau2,
                             beta_hat,
                             X,
                             W,
                             D,
                             eig_vals,
                             model_type) {
  # (-n / 2) log tau^2
  term1 <- (-0.5 * nrow(X)) * log(tau2)

  # eta - X beta
  vector_term <- eta_hat - (X %*% beta_hat)
  # (eta - X beta)^\top Q (eta - X beta)
  term3 <- as.numeric(crossprod(Q %*% vector_term, vector_term))
  # (-1/2 tau^2) x above
  term3 <- (-0.5 / tau2) * term3

  # Trace[Q V]
  # Trace tricks: Tr(A B) = <vec(A), vec(B^T)>
  # For sparse Q this is actually faster than crossprod?
  term4 <- sum(Q * t(Vhat), na.rm = TRUE)
  # term4 <- crossprod(as.vector(Q), as.vector(t(Vhat)))
  # (-1/2 tau^2) x above
  term4 <- (-0.5 / tau2) * term4

  # (1/2) log det Q
  # CAR: Q = D - gamma W = D[I - gamma Z]
  # Log det Q = log det D + log det[I - gamma Z]
  # Given eigenvalues lambda_i of Z, eigenvalues of I - gamma Z are
  # 1 - gamma lambda_i
  # Hence, log det D = sum log D_{i, i}
  # and log det[I - gamma Z] = sum log(1 - gamma lambda_i)
  if ("CAR" == model_type) {
    log_det_D <- sum(log(diag(D)))
    log_det_IgZ <- sum(log(1.0 - gamma_hat * eig_vals), na.rm = TRUE)
    term2 <- 0.5 * (log_det_D + log_det_IgZ)
  }
  # SAR: Q = [I - gamma Z]^\top D [I - gamma Z]
  # log det Q = 2 log det[I - gamma Z] + log det D
  # Similar to CAR, different linear combination
  else if ("SAR" == model_type) {
    log_det_D <- sum(log(diag(D)))
    log_det_IgZ <- sum(log(1.0 - gamma_hat * eig_vals), na.rm = TRUE)
    term2 <- 0.5 * (log_det_D + 2 * log_det_IgZ)
  }
  # Leroux: Q = gamma (D - W) + (1 - gamma) I
  # Given eigenvalues kappa_i of D - W, eigenvalues of Q are
  # gamma kappa_i + (1 - gamma) = gamma (kappa_i - 1) + 1
  # and log det Q = sum log(gamma (kappa_i - 1) + 1)
  else if ("Leroux" == model_type) {
    term2 <- 0.5 * sum(log(gamma_hat * (eig_vals - 1) + 1), na.rm = TRUE)
  }
  # Eigenvalues of Q are passed in.
  else if ("spNNGP" == model_type) {
    term2 <- 0.5 * sum(log(eig_vals), na.rm = TRUE)
  }
  else {
    stop("Invalid model_type")
  }

  return(term1 + term2 + term3 + term4)
}


## Sparse NN GP Kernels

#' Matern Kernel Function.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param d Vector of distances.
#' @param sigma2 Scale parameter.
#' @param rho Range parameter.
#' @param nu (kappa) Order of kernel.
#'
#' @returns A vector of covariance values.
kernel.matern <- function(d, sigma2, rho, nu) {
  # Scaled distance
  z <- (sqrt(2.0 * nu) / rho) * abs(d)

  # sigma^2 2^(1 - nu) / gamma(nu), but logged
  constant_term <- log(sigma2) + (1.0 - nu) * log(2.0) - lgamma(nu)
  # z^nu but logged
  z_term <- nu * log(z)
  # Exponentiate to undo log and multiple with Bessel
  result <- exp(constant_term + z_term) * besselK(z, nu)
  result[d == 0.0] <- sigma2 # Handle zero distances

  return(result)
}

#' Exponential Kernel Function.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param d Vector of distances.
#' @param sigma2 Scale parameter.
#' @param rho Range parameter.
#'
#' @returns A vector of covariance values.
kernel.exp <- function(d, sigma2, rho) {
  return(sigma2 * exp(-abs(d / rho)))
}

#' Gaussian (Squared Exponential) Kernel Function.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param d Vector of distances.
#' @param sigma2 Scale parameter.
#' @param rho Range parameter.
#'
#' @returns A vector of covariance values.
kernel.gauss <- function(d, sigma2, rho) {
  return(sigma2 * exp(-0.5 * (d / rho)^2))
}

#' Spherical Kernel Function.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param d Vector of distances.
#' @param sigma2 Scale parameter.
#' @param rho Range parameter.
#'
#' @returns A vector of covariance values.
kernel.sph <- function(d, sigma2, rho) {
  # Scaled distance values
  z <- abs(d / rho)
  # Cap at 1.0
  z <- pmin(z, 1)

  # 1 - (3/2) z + (1/2) z^3
  result <- 1.0 - 1.5 * z + 0.5 * z^3
  result <- sigma2 * result

  return(result)
}

#' Form a sparse Nearest-Neighbor Gaussian Process Precision Matrix.
#' See https://mc-stan.org/users/documentation/case-studies/nngp.html.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param sp_dist Output of sparseDist_LT function.
#'  Basically, top few rows are nearest distances, bottom few rows are indices.
#'  CALL THAT FUNCTION OR SEE IT FOR MORE DETAILS.
#' @param coords Matrix (x, y) of coordinates.
#' @param cov_type String for covariance model type.
#'  "Exp", "Sph", "Gau", and "Mat" are supported.
#' @param cov_params Covariance/kernel parameters.
#'  Same order as in variogram functions.
#'  Nugget variance (partial sill), Spatial variance (partial sill),
#'  Spatial range (scales distance); if Matern, smoothness kappa.
#'
#' @return A list with a sparse precision matrix and associated eigenvalues.
#' @returns Sparse precision matrix Q.
#' @returns Eigenvalues Dinv.
#' @returns Lower triangular factor A.
#'
#' @note Requires the Matrix library.
#'
#' @import Matrix
#' @importFrom Matrix Diagonal
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix solve
#' @export
nngp_prec_mat <- function(sp_dist, coords, cov_type, cov_params) {
  # Assumed format of sp_dist
  nngp_k <- nrow(sp_dist) / 2

  # Entries of lower triangular matrix A
  A_col <- list()
  A_row <- list()
  A_val <- list()
  # Entries of diagonal D^{-1}
  Dinv_val <- rep(0, ncol(sp_dist) + 1)
  # Loop over samples
  for (idx in 1:ncol(sp_dist)) {
    keep_idx <- which(!is.na(sp_dist[(1 + nngp_k):(2 * nngp_k), idx]))
    if (0 == length(keep_idx)) {
      next
    }

    # C(point, neighbors) vector
    cov_nb <- sp_dist[keep_idx, idx] # Extract distances
    # C(neighbors, neighbors) matrix
    if (1 < length(keep_idx)) {
      cov_nb_mat <- as.matrix(dist(coords[sp_dist[nngp_k + keep_idx, idx], ], diag =
                                     TRUE, upper = TRUE))
    } else {
      cov_nb_mat <- as.matrix(0.0)
    }

    # Call spatial covariance functions
    if ("Mat" == cov_type) {
      cov_nb <- kernel.matern(cov_nb, cov_params[2], cov_params[3], cov_params[4])
      cov_nb_mat <- kernel.matern(cov_nb_mat, cov_params[2], cov_params[3], cov_params[4])
    } else if ("Exp" == cov_type) {
      cov_nb <- kernel.exp(cov_nb, cov_params[2], cov_params[3])
      cov_nb_mat <- kernel.exp(cov_nb_mat, cov_params[2], cov_params[3])
    } else if ("Gau" == cov_type) {
      cov_nb <- kernel.gauss(cov_nb, cov_params[2], cov_params[3])
      cov_nb_mat <- kernel.gauss(cov_nb_mat, cov_params[2], cov_params[3])
    } else if ("Sph" == cov_type) {
      cov_nb <- kernel.sph(cov_nb, cov_params[2], cov_params[3])
      cov_nb_mat <- kernel.sph(cov_nb_mat, cov_params[2], cov_params[3])
    } else {
      stop("Invalid or unimplemented covariance structure.")
    }

    # (C(nb, nb) + tau^2)^{-1} C(point, nb)
    vector_term <- Matrix::solve(cov_nb_mat + diag(cov_params[1], length(cov_nb)), cov_nb)

    # Store elements in A
    A_row[[idx]] <- rep(idx + 1, length(keep_idx))
    A_col[[idx]] <- sp_dist[nngp_k + keep_idx, idx]
    A_val[[idx]] <- as.vector(vector_term)

    # D: C(point, point) + tau^2 - C(point, nb)^\top (C(nb, nb) + tau^2)^{-1} C(point, nb)
    Dinv_val[idx] <- 1.0 / (cov_params[1] + cov_params[2] - crossprod(cov_nb, vector_term))
  }
  # Add last point
  Dinv_val[length(Dinv_val)] <- 1.0 / (cov_params[1] + cov_params[2])

  # Create lower triangular matrix A
  A_col <- unlist(A_col)
  A_row <- unlist(A_row)
  A_val <- unlist(A_val)
  A_mat <- Matrix::sparseMatrix(
    i = A_row,
    j = A_col,
    x = A_val,
    dims = 1 + c(ncol(sp_dist), ncol(sp_dist)),
    triangular = TRUE
  )

  # Identity matrix
  id_mat <- Matrix::Diagonal(dim(A_mat)[1], 1)

  # (I - A)^\top D^{-1} (I - A)
  Q <- crossprod(id_mat - A_mat,
                 Matrix::Diagonal(length(Dinv_val), Dinv_val) %*% (id_mat - A_mat))

  return(list(Q = Q, Dinv = Dinv_val, A = A_mat))
}
