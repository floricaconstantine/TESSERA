## Functions for simulations, e.g., data generation.
# Dependencies in file: Matrix.
# Dependencies: Functions from models.R

#' Generate data for one area.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param n_points Number of points in the area.
#' @param nb_dist Distance threshold for determining if two points are neighbors.
#'  Coordinates are samples in (0, 1)^2 so 0.03 is a reasonable value.
#' @param model_type "CAR", "SAR", or "Leroux"---model for random effects.
#' @param beta_true True fixed effects.
#' @param gamma_true True correlation parameter.
#' @param tau2_true True scale parameter.
#' @param X_type How to generate the covariates X.
#'  Note that the coordinates are sorted by x then by y, so nearby points are
#'  more likely to be neighbors.
#'  "rand_bern": X_{i, j} ~ 2 ( sign(N(0, 1)) + 1 ) in {0, 1}
#'  "rand_unif": X_{i, j} ~ U(0, 1)
#'  "rand_norm": X_{i, j} ~ N(0, 1)
#'  "ar1": X_{:, j} follows an AR(1) process.
#'  "ar1_bern": Generate the AR(1) as above, and apply the same sign transformation
#'    as in "rand_bern"
#'  "intercept": Constant column of all 1s
#' @param ar_gamma Autocorrelation parameter for X if needed.
#'
#' @return A list with the following fields:
#' @returns X: Binary covariates.
#' @returns W: Adjacency matrix.
#' @returns D: Diagonal degree matrix (row-sums of W).
#' @returns eig_list: List of eigenvalues of Z = D^{-1} W (CAR/SAR) or D - W (Leroux).
#' @returns library_size: Scaling for theta.
#' @returns x_coords: x-coordinates of points.
#' @returns y_coords: y-coordinates of points.
#' @returns Q: Unscaled precision matrix. Depends on W, D, gamma_true.
#' @returns phi_true: True sampled random effects.
#'  (multivariate normal with mean zero and covariance tau^2 Q^{-1}).
#' @returns eta_true: phi + X beta.
#' @returns theta_true: exp(eta).
#' @returns z: Sampled counts Pois(theta x lib size).
#'
#' @note Requires the Matrix library.
#' @note Calls the various Q_matrix creation functions.
#'
#' @import Matrix
#' @importFrom Matrix Matrix
#' @importFrom Matrix chol
#' @importFrom Matrix rowSums
#' @importFrom Matrix solve
#' @importFrom Matrix t
#' @importFrom stats arima.sim
#' @importFrom stats dist
#' @importFrom stats rnorm
#' @importFrom stats rpois
#' @importFrom stats runif
#'
#' @export
generate_data_one_area <- function(n_points,
                                   nb_dist,
                                   model_type,
                                   beta_true,
                                   gamma_true = NULL,
                                   tau2_true = NULL,
                                   X_type = "rand_bern",
                                   ar_gamma = 0.75) {
  # Model covariance function
  if ("CAR" == model_type) {
    Q_fcn = Q_matrix_CAR
  } else if ("SAR" == model_type) {
    Q_fcn = Q_matrix_SAR
  } else if ("Leroux" == model_type) {
    Q_fcn = Q_matrix_Leroux
  } else {
    stop("Invalid model_type")
  }

  # Covariate matrix
  if ("rand_bern" == X_type) {
    X <- matrix((sign(stats::rnorm(
      length(beta_true) * n_points
    )) + 1) / 2, nrow = n_points)
  } else if ("rand_norm" == X_type) {
    X <- matrix(stats::rnorm(length(beta_true) * n_points), nrow = n_points)
  } else if ("rand_unif" == X_type) {
    X <- matrix(stats::runif(length(beta_true) * n_points), nrow = n_points)
  } else if (("ar1" == X_type) || ("ar1_bern" == X_type)) {
    X <- matrix(NA, nrow = n_points, ncol = length(beta_true))
    for (col_idx in 1:ncol(X)) {
      X[, col_idx] <- stats::arima.sim(list(ar = ar_gamma), n_points)
    }

    if ("ar1_bern" == X_type) {
      X <- (sign(X) + 1) / 2
    }
  } else if ("intercept" == X_type) {
    X <- matrix(data = 1, nrow = n_points)
  }


  # Expected value/offset for each point
  library_size <- rep(1.0, n_points)

  # Sample x and y coordinates for locations
  x_coords <- stats::runif(n_points)
  y_coords <- stats::runif(n_points)
  coords <- cbind(x_coords, y_coords)
  # Sort by x then y for easier/faster indexing
  coords <- coords[order(coords[, 1], coords[, 2]), ]

  # Compute distance matrix
  d_mat <- as.matrix(stats::dist(coords))
  # Adjacency matrix
  W <- 1 * (d_mat < nb_dist)
  W <- W - diag(diag(W))

  # Assign neighbors to isolated points
  zero_indices <- as.numeric(which(rowSums(W) == 0))
  for (idx in 1:length(zero_indices)) {
    set_idx <- order(d_mat[zero_indices[idx], ])[2]
    W[zero_indices[idx], set_idx] <- 1
    W[set_idx, zero_indices[idx]] <- 1
  }

  # Degree matrix
  D <- diag(rowSums(W))

  # Make sparse
  D <- Matrix::Matrix(D, sparse = TRUE)
  W <- Matrix::Matrix(W, sparse = TRUE)

  # D^{-1} W for CAR and SAR, D - W for Leroux
  if (("CAR" == model_type) || ("SAR" == model_type)) {
    eig_list <-
      Re(eigen(Matrix::solve(D, W), FALSE, only.values = TRUE)$values)
  } else if ("Leroux" == model_type) {
    eig_list <-
      Re(eigen(D - W, TRUE, only.values = TRUE)$values)
  }

  # Check that W is symmetric and has no isolated points
  stopifnot(0 == sum(abs(W - Matrix::t(W))))
  stopifnot(0 == sum(Matrix::rowSums(W) == 0))

  # Generate inverse precision (unscaled) of random effects
  Q <- Q_fcn(W, D, gamma_true)

  # Spatial random effects
  phi_true <- as.numeric(sqrt(tau2_true)
                         * Matrix::solve(Matrix::chol(Q), stats::rnorm(n_points)))
  # Add in covariate effect
  eta_true <- as.numeric(X %*% beta_true) + phi_true
  # Get Poisson parameter
  theta_true <- exp(eta_true)
  # Generate data
  z <- stats::rpois(n_points, theta_true * library_size)

  return(
    list(
      X = X,
      W = W,
      D = D,
      eig_list = eig_list,
      library_size = library_size,
      x_coords = x_coords,
      y_coords = y_coords,
      Q = Q,
      phi_true = phi_true,
      eta_true = eta_true,
      theta_true = theta_true,
      z = z
    )
  )
}


#' Generate data for one area from the spNNGP model.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param n_points Number of points in the area.
#' @param nb_dist Distance threshold for determining if two points are neighbors.
#'  Coordinates are samples in (0, 1)^2 so 0.03 is a reasonable value.
#' @param cov_type Which model for the Gaussian process covariance.
#'    "Exp", "Mat", "Gau", and "Sph" are the valid options, for
#'    Exponential, Matern, Gaussian, and Spherical, respectively.
#' @param cov_params Variogram/covariance parameters.
#' @param nngp_k Number of nearest neighbors.
#' @param beta_true True fixed effects.
#' @param X_type How to generate the covariates X.
#'  Note that the coordinates are sorted by x then by y, so nearby points are
#'  more likely to be neighbors.
#'  "rand_bern": X_{i, j} ~ 2 ( sign(N(0, 1)) + 1 ) in {0, 1}
#'  "rand_unif": X_{i, j} ~ U(0, 1)
#'  "rand_norm": X_{i, j} ~ N(0, 1)
#'  "ar1": X_{:, j} follows an AR(1) process.
#'  "ar1_bern": Generate the AR(1) as above, and apply the same sign transformation
#'    as in "rand_bern"
#'  "intercept": Constant column of all 1s
#' @param ar_gamma Autocorrelation parameter for X if needed.
#'
#' @return A list with the following fields:
#' @returns X: Binary covariates.
#' @returns W: Adjacency matrix.
#' @returns D: Diagonal degree matrix (row-sums of W).
#' @returns eig_list: List of eigenvalues of Z = D^{-1} W (CAR/SAR) or D - W (Leroux).
#' @returns library_size: Scaling for theta.
#' @returns x_coords: x-coordinates of points.
#' @returns y_coords: y-coordinates of points.
#' @returns Sparse precision matrix Q.
#' @returns Eigenvalues Dinv.
#' @returns Lower triangular factor A.
#' @returns phi_true: True sampled random effects.
#'  (multivariate normal with mean zero and covariance tau^2 Q^{-1}).
#' @returns eta_true: phi + X beta.
#' @returns theta_true: exp(eta).
#' @returns z: Sampled counts Pois(theta x lib size).
#'
#' @note Requires the Matrix library.
#' @note Calls the nngp_prec_mat creation functions.
#'
#' @import Matrix
#' @importFrom Matrix Matrix
#' @importFrom Matrix chol
#' @importFrom Matrix Diagonal
#' @importFrom Matrix rowSums
#' @importFrom Matrix solve
#' @importFrom Matrix t
#' @importFrom stats arima.sim
#' @importFrom stats dist
#' @importFrom stats rnorm
#' @importFrom stats rpois
#' @importFrom stats runif
#'
#' @export
generate_data_one_area_spNNGP <- function(n_points,
                                          nb_dist,
                                          cov_type,
                                          cov_params,
                                          nngp_k,
                                          beta_true,
                                          X_type = "rand_bern",
                                          ar_gamma = 0.75) {
  # Covariate matrix
  if ("rand_bern" == X_type) {
    X <- matrix((sign(stats::rnorm(
      length(beta_true) * n_points
    )) + 1) / 2, nrow = n_points)
  } else if ("rand_norm" == X_type) {
    X <- matrix(stats::rnorm(length(beta_true) * n_points), nrow = n_points)
  } else if ("rand_unif" == X_type) {
    X <- matrix(stats::runif(length(beta_true) * n_points), nrow = n_points)
  } else if (("ar1" == X_type) || ("ar1_bern" == X_type)) {
    X <- matrix(NA, nrow = n_points, ncol = length(beta_true))
    for (col_idx in 1:ncol(X)) {
      X[, col_idx] <- stats::arima.sim(list(ar = ar_gamma), n_points)
    }

    if ("ar1_bern" == X_type) {
      X <- (sign(X) + 1) / 2
    }
  } else if ("intercept" == X_type) {
    X <- matrix(data = 1, nrow = n_points)
  }

  # Expected value/offset for each point
  library_size <- rep(1.0, n_points)

  # Sample x and y coordinates for locations
  x_coords <- stats::runif(n_points)
  y_coords <- stats::runif(n_points)
  coords <- cbind(x_coords, y_coords)
  # Sort by x then y for easier/faster indexing
  coords <- coords[order(coords[, 1], coords[, 2]), ]

  # Compute distance matrix
  d_mat <- as.matrix(stats::dist(coords))
  # Adjacency matrix
  W <- 1 * (d_mat < nb_dist)
  W <- W - diag(diag(W))

  # Assign neighbors to isolated points
  zero_indices <- as.numeric(which(rowSums(W) == 0))
  for (idx in 1:length(zero_indices)) {
    set_idx <- order(d_mat[zero_indices[idx], ])[2]
    W[zero_indices[idx], set_idx] <- 1
    W[set_idx, zero_indices[idx]] <- 1
  }

  # Degree matrix
  D <- diag(rowSums(W))

  # Make sparse
  D <- Matrix::Matrix(D, sparse = TRUE)
  W <- Matrix::Matrix(W, sparse = TRUE)

  # Check that W is symmetric and has no isolated points
  stopifnot(0 == sum(abs(W - Matrix::t(W))))
  stopifnot(0 == sum(Matrix::rowSums(W) == 0))

  # Generate inverse precision (unscaled) of random effects
  sp_dist <- sparseDist_LT(coords, nngp_k)
  close_to_zero_const <- 2.0 * .Machine$double.eps
  if ((close_to_zero_const >= cov_params[1]) &&
      (close_to_zero_const >= cov_params[2])) {
    print("BOTH NUGGET AND SPATIAL VARIANCE ARE ZERO.")
    print("FAILSAFE: Q IS ZERO AND NOT USED")

    Q <- Matrix::Diagonal(n_points, 0)
    A <- Matrix::Diagonal(n_points, 0)
    Dinv <- Matrix::Diagonal(n_points, 0)
  } else {
    param_est <- nngp_prec_mat(sp_dist, coords, cov_type, cov_params)
    Q <- param_est$Q
    A <- param_est$A
    Dinv <- param_est$Dinv
  }

  # Spatial random effects
  close_to_zero_const <- 2.0 * .Machine$double.eps
  if ((close_to_zero_const >= cov_params[1]) &&
      (close_to_zero_const >= cov_params[2])) {
    print("BOTH NUGGET AND SPATIAL VARIANCE ARE ZERO.")
    print("FAILSAFE: PHI IS ZERO.")
    eta_true <- as.numeric(X %*% beta_true)
    phi_true <- 0 * eta_true
  } else{
    phi_true <- as.numeric(Matrix::solve(Matrix::chol(Q), stats::rnorm(n_points)))
    # Add in covariate effect
    eta_true <- as.numeric(X %*% beta_true) + phi_true
  }

  # Get Poisson parameter
  theta_true <- exp(eta_true)
  # Generate data
  z <- stats::rpois(n_points, theta_true * library_size)

  return(
    list(
      X = X,
      W = W,
      D = D,
      library_size = library_size,
      x_coords = x_coords,
      y_coords = y_coords,
      Q = Q,
      A = A,
      Dinv = Dinv,
      phi_true = phi_true,
      eta_true = eta_true,
      theta_true = theta_true,
      z = z
    )
  )
}
