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
#' @param precomp Optional precomputed matrices (unused for CAR).
#'
#' @returns Unscaled precision matrix.
Q_matrix_CAR <- function(W, D, gamma_val, precomp = NULL) {
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
#' @param precomp Precomputed W D^\{-1\} W matrix.
#'
#' @returns Unscaled precision matrix.
#'
#' @note Requires the Matrix library.
#'
#' @import Matrix
#' @importFrom Matrix Diagonal
#' @importFrom Matrix diag
Q_matrix_SAR <- function(W, D, gamma_val, precomp = NULL) {
  if (is.null(precomp)) {
    D_inv <- Matrix::Diagonal(n = dim(W)[1], x = 1 / Matrix::diag(D))
    precomp <- W %*% D_inv %*% W
  }
  # Q = D - 2*gamma*W + gamma^2 * W D^{-1} W
  return(D - (2 * gamma_val * W) + ((gamma_val^2) * precomp))
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
#' @param precomp Precomputed list containing D_minus_W and id_mat.
#'
#' @returns Unscaled precision matrix.
#'
#' @note Requires the Matrix library.
#'
#' @import Matrix
Q_matrix_Leroux <- function(W, D, gamma_val, precomp = NULL) {
  if (is.null(precomp)) {
    precomp <- list(D_minus_W = D - W,
                    id_mat = Matrix::Diagonal(n = dim(W)[1], 1))
  }
  # gamma (D - W) + (1 - gamma) I
  return((gamma_val * precomp$D_minus_W) + ((1 - gamma_val) * precomp$id_mat))
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
#'
#' @importFrom Matrix diag
#' @importFrom methods as
neg_hessian_tau2 <- function(Vhat, eta_hat, Q, tau2_hat, beta_hat, X) {
  # eta - X beta
  vector_term <- eta_hat - (X %*% beta_hat)
  
  # term1: crossprod(Q %*% vector_term, vector_term) is already efficient
  term1 <- as.numeric(crossprod(Q %*% vector_term, vector_term))
  
  # term2: Trace[Q V]
  # We extract the row/column indices of the non-zero elements of Q
  # This avoids ever looking at the dense off-diagonal elements of Vhat
  # Old: Q_sparse <- methods::as(Q, "dgCMatrix")
  Q_sparse <- methods::as(methods::as(Q, "generalMatrix"), "CsparseMatrix")
  idx_Q <- cbind(Q_sparse@i + 1L, rep(1:ncol(Q_sparse), diff(Q_sparse@p)))
  
  # We only compute Vhat[i, j] for non-zero Q[i, j]
  # Vhat is symmetric, so Vhat[i, j] == Vhat[j, i]
  term2 <- sum(Q_sparse@x * Vhat[idx_Q])
  
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
#' @param eig_vals Eigenvalues of Z = D^\{-1\} W.
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
#' @param eig_vals Eigenvalues of Z = D^\{-1\} W.
#'
#' @returns Negative Hessian of gamma.
#'  Also called the observed information.
#'
#' @note Requires the Matrix library.
#'
#' @import Matrix
#' @importFrom Matrix Diagonal
#' @importFrom Matrix diag
#' @importFrom methods as
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
  zeta2 <- as.numeric(crossprod(W %*% vector_term, Z %*% vector_term))
  
  # Tr[W Z V]
  # Create sparse matrix WZ to avoid dense matrix allocations
  # Old: WZ_sparse <- methods::as(W %*% Z, "dgCMatrix")
  WZ_sparse <- methods::as(methods::as(W %*% Z, "generalMatrix"), "CsparseMatrix")
  
  # Extract row/col indices of non-zero elements
  idx_WZ <- cbind(WZ_sparse@i + 1L, rep(1:ncol(WZ_sparse), diff(WZ_sparse@p)))
  
  # Sum the element-wise product. Since Vhat is symmetric, Vhat[i, j] == Vhat[j, i]
  zeta4 <- sum(WZ_sparse@x * Vhat[idx_WZ])
  
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
#' @param Vhat Estimated covariance matrix of eta (sparse subset).
#' @param eta_hat Estimated mean of eta.
#' @param Q Unscaled precision matrix.
#'  The scaled precision matrix in spNNGP.
#' @param gamma_hat Parameter of interest: Estimate of correlation parameter.
#' @param tau2 Parameter of interest: Current estimate of scale.
#' @param beta_hat Current estimate of covariate effects.
#' @param X Covariate matrix.
#' @param W Neighbor/adjacency matrix (symmetric, binary).
#' @param D Degree matrix (diagonal, values are row-sums of W).
#' @param eig_vals Eigenvalues of Z = D^\{-1\} W (CAR/SAR), D - W (Leroux), or Q (spNNGP).
#' @param model_type "CAR", "SAR", or "Leroux", or "spNNGP".
#'    If spNNGP, set tau^2 to 1 and the value of gamma is irrelevant.
#'  Model for the random effects.
#'
#' @returns Expected log likelihood.
#'
#' @importFrom methods as
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
  # Convert both Q and Vhat to general CSC (dgCMatrix) and compute the trace
  # in C++ via a linear two-pointer scan over their CSC slots
  Q_sparse <- methods::as(methods::as(Q, "generalMatrix"), "CsparseMatrix")
  V_sparse <- methods::as(methods::as(Vhat, "generalMatrix"), "CsparseMatrix")
  tr_QV <- fast_trace_quad_cpp(Q_sparse@i,
                               Q_sparse@p,
                               Q_sparse@x,
                               V_sparse@i,
                               V_sparse@p,
                               V_sparse@x)
  term4 <- (-0.5 / tau2) * tr_QV
  
  # (1/2) log det Q
  if ("CAR" == model_type) {
    log_det_D <- sum(log(Matrix::diag(D)))
    log_det_IgZ <- sum(log(1.0 - gamma_hat * eig_vals), na.rm = TRUE)
    term2 <- 0.5 * (log_det_D + log_det_IgZ)
  }
  else if ("SAR" == model_type) {
    log_det_D <- sum(log(Matrix::diag(D)))
    log_det_IgZ <- sum(log(1.0 - gamma_hat * eig_vals), na.rm = TRUE)
    term2 <- 0.5 * (log_det_D + 2 * log_det_IgZ)
  }
  else if ("Leroux" == model_type) {
    term2 <- 0.5 * sum(log(gamma_hat * (eig_vals - 1) + 1), na.rm = TRUE)
  }
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
kernel_matern <- function(d, sigma2, rho, nu) {
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
kernel_exp <- function(d, sigma2, rho) {
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
kernel_gauss <- function(d, sigma2, rho) {
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
kernel_sph <- function(d, sigma2, rho) {
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
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param sp_dist Output of sparse_dist_LT function.
#' @param nb_dist Precomputed list of Euclidean distance matrices for neighbors.
#' @param cov_type String for covariance model type.
#' @param cov_params Covariance/kernel parameters. Assuming gstat conventions,
#'  nugget (non-spatial), sill (spatial/scale), range, order (Matern).
#'
#' @returns A list with a sparse precision matrix Q, eigenvalues Dinv, and factor A.
#'
#' @import Matrix
nngp_prec_mat <- function(sp_dist, nb_dist, cov_type, cov_params) {
  n_cols <- ncol(sp_dist)
  
  # Compute triplet vectors (A_row, A_col, A_val) and Dinv in C++ via LAPACK dgesv
  cpp_out <- nngp_prec_triplets_cpp(sp_dist, nb_dist, cov_type, cov_params)
  
  A_mat <- Matrix::sparseMatrix(
    i = cpp_out$A_row,
    j = cpp_out$A_col,
    x = cpp_out$A_val,
    dims = 1 + c(n_cols, n_cols),
    triangular = TRUE
  )
  
  Dinv_val <- cpp_out$Dinv
  id_mat <- Matrix::Diagonal(dim(A_mat)[1], 1)
  I_minus_A <- id_mat - A_mat
  Q <- Matrix::crossprod(I_minus_A, Matrix::Diagonal(length(Dinv_val), Dinv_val) %*% I_minus_A)
  
  return(list(Q = Q, Dinv = Dinv_val, A = A_mat))
}
