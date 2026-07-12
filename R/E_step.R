## E-Step (Expectation Step) functions.
# Dependencies in file: Matrix, pracma.

#' Compute the Covariance Matrix of the random effects eta.
#' Part of the E-Step in the EM algorithm.
#' Follows Clayton and Kaldor (1987) for the Poisson-Gaussian model.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @note This method is an Empirical Bayes approximation.
#' @note This method inverts a matrix, because the inverse is actually needed...
#'
#' @param Q Unscaled precision matrix.
#' @param tau2 Precision/covariance matrix scaling.
#' @param z Observed counts.
#'
#' @returns Estimated covariance matrix.
#'
#' @note Requires the Matrix library.
#' @note Requires the pracma library.
#'
#' @import Matrix
#' @importFrom Matrix solve
#' @importFrom Matrix Diagonal
#' @importFrom pracma pinv
E_step_Vhat <- function(Q, tau2, z) {
  # V^{-1} = Q / tau^2 + Diagonal(0.5 + z)
  Vinv <- (Q / tau2) + Matrix::Diagonal(dim(Q)[1], 0.5 + z)
  
  # Compute V
  # Try basic inversion first
  err_flag <- tryCatch({
    V <- Matrix::solve(Vinv)
    return(V)
  }, error = function(cond) {
    warning("INVERSION OF Inv(V_hat) FAILED")
    err_flag <- TRUE
  })
  # Then try pseudoinverse
  if (err_flag) {
    err_flag <- tryCatch({
      V <- pracma::pinv(as.matrix(Vinv))
      return(V)
    }, error = function(cond) {
      warning("PSEUDOINVERSE OF Inv(V_hat) FAILED")
      err_flag <- TRUE
    })
  }
  if (err_flag) {
    stop("PSEUDOINVERSE OF Inv(V_hat FAILED: CANNOT CONTINUE")
  }
}

#' Compute the INVERSE Covariance Matrix of the random effects eta.
#' Part of the E-Step in the EM algorithm.
#' Follows Clayton and Kaldor (1987) for the Poisson-Gaussian model.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @note This method is an Empirical Bayes approximation.
#'
#' @param Q Unscaled precision matrix.
#' @param tau2 Precision/covariance matrix scaling.
#' @param z Observed counts.
#'
#' @returns Estimated covariance matrix P LU decomposition.
#'
#' @note Requires the Matrix library.
#'
#' @import Matrix
#' @importFrom Matrix Diagonal
#' @importFrom Matrix lu
#' @importFrom Matrix expand
E_step_Vhat_PLU <- function(Q, tau2, z) {
  # V^{-1} = Q / tau^2 + Diagonal(0.5 + z)
  Vinv <- (Q / tau2) + Matrix::Diagonal(dim(Q)[1], 0.5 + z)
  
  # Compute PLU decomposition: V^{-1} = P L U
  decomp  <- Matrix::expand(Matrix::lu(Vinv))
  decomp$Linv <- Matrix::solve(decomp$L)
  decomp$Uinv <- Matrix::solve(decomp$U)
  return(decomp)
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
#' @param Vhat Estimated covariance matrix of eta.
#' @param Q Unscaled precision matrix.
#' @param tau2 Precision/covariance matrix scaling.
#' @param beta_hat Current estimate of covariate effects.
#' @param X Covariate matrix.
#' @param z Observed counts.
#' @param N Library size or expected counts.
#'
#' @returns Estimated mean of eta.
E_step_etahat <- function(Vhat, Q, tau2, beta_hat, X, z, N) {
  # (z_i + 1/2) log[(z_i + 1/2) / N_i] - (1/2)
  term1 <- z + 0.5
  term1 <- term1 * log(term1 / N) - 0.5
  
  # (Q / tau^2) X beta
  term2 <- (Q %*% (X %*% beta_hat)) / tau2
  
  # V [term1 + term2]
  eta_hat <- Vhat %*% (term1 + term2)
  
  return(eta_hat)
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
#' @param decomp Estimated covariance matrix of eta PLU decomposition.
#' @param Q Unscaled precision matrix.
#' @param tau2 Precision/covariance matrix scaling.
#' @param beta_hat Current estimate of covariate effects.
#' @param X Covariate matrix.
#' @param z Observed counts.
#' @param N Library size or expected counts.
#'
#' @returns Estimated mean of eta.
#'
#' @importFrom Matrix solve
E_step_etahat_PLU <- function(decomp, Q, tau2, beta_hat, X, z, N) {
  # (z_i + 1/2) log[(z_i + 1/2) / N_i] - (1/2)
  term1 <- z + 0.5
  term1 <- term1 * log(term1 / N) - 0.5
  
  # (Q / tau^2) X beta
  term2 <- (Q %*% (X %*% beta_hat)) / tau2
  
  # V [term1 + term2]
  # eta_hat <- Matrix::solve(Vhat_inv, (term1 + term2))
  eta_hat <- decomp$Uinv %*% (decomp$Linv %*% (t(decomp$P) %*% (term1 + term2)))
  
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
#' @param Vhat Estimated covariance matrix of eta.
#' @param eta_hat Estimated mean of eta.
#'
#' @returns Estimated poisson parameters theta_hat.
E_step_thetahat <- function(Vhat, eta_hat) {
  # exp(eta + (1/2) V_{i, i})
  return(exp(eta_hat + 0.5 * diag(Vhat)))
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
#' @param decomp Estimated covariance matrix of eta PLU decomposition.
#' @param eta_hat Estimated mean of eta.
#'
#' @returns Estimated poisson parameters theta_hat.
E_step_thetahat_PLU <- function(decomp, eta_hat) {
  # exp(eta + (1/2) V_{i, i})
  # Take an P L U decomposition and use that to get the diagonals.
  # Inverse of P L U = U^{-1} L^{-1} P^T
  diag_lu <- colSums((decomp$Linv %*% t(decomp$P)) * t(decomp$Uinv))
  
  return(exp(eta_hat + 0.5 * diag_lu))
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
