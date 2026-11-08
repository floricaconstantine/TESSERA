## Functions specific to the models/distributions used.
# Dependencies in file: Matrix.

## Covariance matrices

#' Compute unscaled precision matrix for CAR model.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param W Neighbor/adjacency matrix (symmetric, binary).
#' @param D Degree matrix (diagonal, values are row-sums of W).
#' @param gamma_hat Spatial correlation parameter.
#' @param precomp Optional precomputed list from \code{prep_Q_lattice_cpp}.
#'
#' @returns Unscaled precision matrix Q.
Q_matrix_CAR <- function(W, D, gamma_hat, precomp = NULL) {
  if (!is.null(precomp) &&
      is.list(precomp) && !is.null(precomp$tmpl)) {
    return(eval_Q_lattice_cpp(precomp, c(1.0, -gamma_hat)))
  }
  # Q = D - gamma W
  return(D - gamma_hat * W)
}

#' Compute unscaled precision matrix for SAR model.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param W Neighbor/adjacency matrix (symmetric, binary).
#' @param D Degree matrix (diagonal, values are row-sums of W).
#' @param gamma_hat Spatial correlation parameter.
#' @param precomp Optional precomputed list from \code{prep_Q_lattice_cpp}
#'  or precomputed matrix \code{W \%*\% D^{-1} \%*\% W}.
#'
#' @returns Unscaled precision matrix Q.
#'
#' @note Requires the Matrix library.
#'
#' @import Matrix
#' @importFrom Matrix Diagonal
#' @importFrom Matrix diag
Q_matrix_SAR <- function(W, D, gamma_hat, precomp = NULL) {
  if (!is.null(precomp) &&
      is.list(precomp) && !is.null(precomp$tmpl)) {
    return(eval_Q_lattice_cpp(precomp, c(1.0, -2.0 * gamma_hat, gamma_hat^2)))
  }
  if (is.null(precomp)) {
    # Fallback if called outside the main EM loop
    D_inv <- Matrix::Diagonal(dim(W)[1], 1 / Matrix::diag(D))
    WZW <- W %*% (D_inv %*% W)
  } else {
    WZW <- precomp
  }
  
  # Q = D - 2 gamma W + gamma^2 W D^{-1} W
  return(D - (2.0 * gamma_hat) * W + (gamma_hat^2) * WZW)
}

#' Compute unscaled precision matrix for Leroux model.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param W Neighbor/adjacency matrix (symmetric, binary).
#' @param D Degree matrix (diagonal, values are row-sums of W).
#' @param gamma_hat Spatial correlation parameter.
#' @param precomp Optional precomputed list from \code{prep_Q_lattice_cpp}
#'  or list containing \code{id_mat} and \code{D_minus_W}.
#'
#' @returns Unscaled precision matrix Q.
#'
#' @note Requires the Matrix library.
#'
#' @import Matrix
#' @importFrom Matrix Diagonal
Q_matrix_Leroux <- function(W, D, gamma_hat, precomp = NULL) {
  if (!is.null(precomp) &&
      is.list(precomp) && !is.null(precomp$tmpl)) {
    return(eval_Q_lattice_cpp(precomp, c(1.0 - gamma_hat, gamma_hat)))
  }
  if (is.null(precomp)) {
    # Fallback if called outside the main EM loop
    id_mat <- Matrix::Diagonal(dim(W)[1], 1)
    D_minus_W <- D - W
    return((1.0 - gamma_hat) * id_mat + gamma_hat * D_minus_W)
  } else {
    # Q = (1 - gamma) I + gamma (D - W)
    return((1.0 - gamma_hat) * precomp$id_mat + gamma_hat * precomp$D_minus_W)
  }
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
  
  # term2: Trace[Q V] via C++ CSC two-pointer scan
  Q_sparse <- if (inherits(Q, "dgCMatrix"))
    Q
  else
    methods::as(methods::as(Q, "generalMatrix"), "CsparseMatrix")
  V_sparse <- if (inherits(Vhat, "dgCMatrix"))
    Vhat
  else
    methods::as(methods::as(Vhat, "generalMatrix"), "CsparseMatrix")
  term2 <- fast_trace_quad_cpp(Q_sparse@i,
                               Q_sparse@p,
                               Q_sparse@x,
                               V_sparse@i,
                               V_sparse@p,
                               V_sparse@x)
  
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
  WZ <- W %*% Z
  WZ_sparse <- if (inherits(WZ, "dgCMatrix"))
    WZ
  else
    methods::as(methods::as(WZ, "generalMatrix"), "CsparseMatrix")
  V_sparse  <- if (inherits(Vhat, "dgCMatrix"))
    Vhat
  else
    methods::as(methods::as(Vhat, "generalMatrix"), "CsparseMatrix")
  zeta4 <- fast_trace_quad_cpp(WZ_sparse@i,
                               WZ_sparse@p,
                               WZ_sparse@x,
                               V_sparse@i,
                               V_sparse@p,
                               V_sparse@x)
  
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
  # Avoid S4 .M2C coercion if Q and Vhat are already dgCMatrix
  Q_sparse <- if (inherits(Q, "dgCMatrix"))
    Q
  else
    methods::as(methods::as(Q, "generalMatrix"), "CsparseMatrix")
  V_sparse <- if (inherits(Vhat, "dgCMatrix"))
    Vhat
  else
    methods::as(methods::as(Vhat, "generalMatrix"), "CsparseMatrix")
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
#' @param build_Q Logical: If \code{TRUE} (default), constructs the sparse
#'  precision matrix \code{Q = (I - A)^T D^{-1} (I - A)} and sparse factor \code{A}.
#'  If \code{FALSE}, skips building \code{Q} and \code{A} and instead computes
#'  \code{XtQeta} and \code{XtQX} directly in C++ from \code{X} and \code{eta}.
#' @param X Optional covariate matrix (N x p) used when \code{build_Q = FALSE}.
#' @param eta Optional random effect vector (length N) used when \code{build_Q = FALSE}.
#'
#' @returns A list with a sparse precision matrix Q, eigenvalues Dinv, and factor A
#'  (plus \code{XtQeta} and \code{XtQX} when \code{build_Q = FALSE}).
#'
#' @import Matrix
nngp_prec_mat <- function(sp_dist,
                          nb_dist,
                          cov_type,
                          cov_params,
                          build_Q = TRUE,
                          X = NULL,
                          eta = NULL) {
  n_cols <- ncol(sp_dist)
  
  # Compute triplet vectors (A_row, A_col, A_val) and Dinv in C++
  cpp_out <- nngp_prec_triplets_cpp(sp_dist, nb_dist, cov_type, cov_params)
  Dinv_val <- cpp_out$Dinv
  
  # Fast path for inner CM iterations (opt_idx < opt_iters):
  # Compute X^T Q eta and X^T Q X directly from (I - A) without forming sparse Q
  if (!build_Q && !is.null(X) && !is.null(eta)) {
    beta_terms <- nngp_beta_terms_cpp(cpp_out$A_row,
                                      cpp_out$A_col,
                                      cpp_out$A_val,
                                      Dinv_val,
                                      X,
                                      eta)
    return(
      list(
        Q = NULL,
        Dinv = Dinv_val,
        A = NULL,
        XtQeta = beta_terms$XtQeta,
        XtQX = beta_terms$XtQX
      )
    )
  }
  
  A_mat <- Matrix::sparseMatrix(
    i = cpp_out$A_row,
    j = cpp_out$A_col,
    x = cpp_out$A_val,
    dims = 1 + c(n_cols, n_cols),
    triangular = TRUE
  )
  
  id_mat <- Matrix::Diagonal(dim(A_mat)[1], 1)
  I_minus_A <- id_mat - A_mat
  Q <- Matrix::crossprod(I_minus_A, Matrix::Diagonal(length(Dinv_val), Dinv_val) %*% I_minus_A)
  
  return(list(Q = Q, Dinv = Dinv_val, A = A_mat))
}
