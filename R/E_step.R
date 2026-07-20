## E-Step (Expectation Step) functions.
# Dependencies in file: Matrix, pracma, sparseinv.

#' Compute the Covariance Matrix of the random effects eta.
#' Part of the E-Step in the EM algorithm.
#' Follows Clayton and Kaldor (1987) for the Poisson-Gaussian model.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @note This method is an Empirical Bayes approximation.
#' @note This method computes a sparse inverse subset using the Takahashi equations
#'  to avoid an O(N^3) dense matrix inversion. It only computes the diagonal and
#'  non-zero off-diagonal elements that match the sparsity pattern of the precision matrix.
#'
#' @param Vinv Precomputed Precision matrix.
#' @param P Cholesky Permutation (precomputed).
#'
#' @returns Estimated covariance matrix (sparse subset).
#'
#' @note Requires the Matrix library.
#' @note Requires the sparseinv library.
#'
#' @import Matrix
#' @importFrom Matrix Diagonal
#' @importFrom Matrix solve
#' @importFrom methods as
#' @importFrom sparseinv Takahashi_Davis
E_step_Vhat <- function(Vinv, P = NULL) {
  # V^{-1} = Q / tau^2 + Diagonal(0.5 + z)
  # Vinv <- (Q / tau2) + Matrix::Diagonal(dim(Q)[1], 0.5 + z)
  
  # Try Takahashi sparse inverse subset first
  res <- tryCatch({
    # Vinv_sparse <- methods::as(methods::as(Vinv, "generalMatrix"), "CsparseMatrix")
    Vinv <- methods::as(Vinv, "dgCMatrix")
    V <- sparseinv::Takahashi_Davis(Vinv, P = P)
    methods::as(V, "symmetricMatrix")
  }, error = function(cond) {
    warning("TAKAHASHI SPARSE INVERSE FAILED; FALLING BACK TO DENSE INVERSION")
    NULL
  })
  if (!is.null(res))
    return(res)
  
  # Fallback 1: Standard dense inversion (O(N^3)) if the sparse package is missing
  res <- tryCatch({
    Matrix::solve(Vinv)
  }, error = function(cond) {
    warning("DENSE INVERSION OF Inv(V_hat) FAILED")
    NULL
  })
  if (!is.null(res))
    return(res)
  
  # Fallback 2: Pseudoinverse (Warning: extremely slow for large N)
  res <- tryCatch({
    pracma::pinv(as.matrix(Vinv))
  }, error = function(cond) {
    warning("PSEUDOINVERSE OF Inv(V_hat) FAILED")
    NULL
  })
  if (!is.null(res))
    return(res)
  
  stop("ALL INVERSIONS OF Inv(V_hat) FAILED: CANNOT CONTINUE")
}

#' Compute the Expectation of the random effects eta.
#' Part of the E-Step in the EM algorithm.
#' Follows Clayton and Kaldor (1987) for the Poisson-Gaussian model.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @note This method is an Empirical Bayes approximation.
#'
#' @param Vinv Precomputed Precision matrix.
#' @param Q Unscaled precision matrix.
#' @param tau2 Precision/covariance matrix scaling.
#' @param beta_hat Current estimate of covariate effects.
#' @param X Covariate matrix.
#' @param z Observed counts.
#' @param N Library size or expected counts.
#'
#' @returns Estimated mean of eta.
#'
#' @importFrom Matrix Diagonal
#' @importFrom Matrix solve
E_step_etahat <- function(Vinv, Q, tau2, beta_hat, X, z, N) {
  # (z_i + 1/2) log[(z_i + 1/2) / N_i] - (1/2)
  term1 <- z + 0.5
  term1 <- term1 * log(term1 / N) - 0.5
  
  # (Q / tau^2) X beta
  term2 <- (Q %*% (X %*% beta_hat)) / tau2
  
  # Reconstruct the sparse precision matrix V^{-1}
  # Vinv <- (Q / tau2) + Matrix::Diagonal(dim(Q)[1], 0.5 + z)
  
  # Solve the sparse linear system V^{-1} eta_hat = (term1 + term2)
  eta_hat <- as.numeric(Matrix::solve(Vinv, (term1 + term2)))
  
  return(eta_hat)
}


#' Compute the estimated poisson parameters theta.
#' Note that eta = log theta.
#' Associated with the E-Step in the EM algorithm, but not needed to run.
#' I.e., this is a utility function.
#' Follows Clayton and Kaldor (1987) for the Poisson-Gaussian model.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param Vhat Estimated covariance matrix of eta (sparse subset).
#' @param eta_hat Estimated mean of eta.
#'
#' @returns Estimated poisson parameters theta_hat.
#'
#' @importFrom Matrix diag
E_step_thetahat <- function(Vhat, eta_hat) {
  # exp(eta + (1/2) V_{i, i})
  return(exp(eta_hat + 0.5 * Matrix::diag(Vhat)))
}

#' Compute the best fit to the data.
#' Associated with the E-Step in the EM algorithm, but not needed to run.
#' I.e., this is a utility function.
#' Follows Clayton and Kaldor (1987) for the Poisson-Gaussian model.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param theta_hat Estimated poisson parameters.
#' @param N Library size or expected counts.
#'
#' @returns Estimated values z (model fits).
E_step_predict <- function(theta_hat, N) {
  # theta N - 1/2
  return(theta_hat * N - 0.5)
}
