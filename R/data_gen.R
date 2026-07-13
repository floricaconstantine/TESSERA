## Functions for simulations, e.g., data generation.
# Dependencies in file: Matrix.
# Dependencies: Functions from models.R

#' Simulate spatial count data for a single area
#'
#' Generates synthetic spatial transcriptomics data based on a Poisson
#' Generalized Linear Mixed Model (\eqn{GLMM}). The function simulates spatial random
#' effects \eqn{\phi} using a \eqn{CAR}, \eqn{SAR}, or \eqn{Leroux} model, generates
#' covariates \eqn{X} based on several stochastic processes, and samples counts
#' \eqn{z} from a Poisson distribution.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param n_points Integer: Number of observations/points in the area.
#' @param nb_dist Numeric: Distance threshold for determining adjacency.
#'   Since coordinates are sampled in \eqn{(0, 1)^2}, a value of 0.03 is
#'   typically appropriate for \eqn{n \approx 1000}.
#' @param model_type Character: The spatial model for random effects.
#'   Options are "CAR", "SAR", or "Leroux".
#' @param beta_true Numeric vector: True fixed effect coefficients \eqn{\beta}.
#' @param gamma_true Numeric: True spatial correlation parameter \eqn{\gamma}.
#' @param tau2_true Numeric: True spatial scale parameter \eqn{\tau^2}.
#' @param X_type Character: Method for generating the design matrix \eqn{X}.
#'   Options include:
#'   \itemize{
#'     \item \strong{"rand_bern"}: Bernoulli covariates \eqn{X_{i, j} \in \{0, 1\}}.
#'     \item \strong{"rand_unif"}: Uniformly distributed \eqn{X_{i, j} \sim U(0, 1)}.
#'     \item \strong{"rand_norm"}: Normally distributed \eqn{X_{i, j} \sim N(0, 1)}.
#'     \item \strong{"ar1"}: \eqn{X} follows a spatial AR(1) process across sorted coordinates.
#'     \item \strong{"ar1_bern"}: AR(1) process transformed into binary values.
#'     \item \strong{"intercept"}: A constant column of 1s.
#'   }
#' @param ar_gamma Numeric: Autocorrelation parameter for \eqn{X} (used if \code{X_type}
#'   is "ar1" or "ar1_bern").
#' @param X Optional matrix: A pre-defined design matrix to use instead of generating one.
#' @param library_size Optional numeric vector: Pre-specified library sizes
#'   (offsets). Defaults to a vector of 1s if not provided.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \strong{X}: The generated design matrix of covariates.
#'   \item \strong{W}: The sparse adjacency matrix.
#'   \item \strong{D}: The diagonal degree matrix.
#'   \item \strong{eig_list}: Eigenvalues of the standardized adjacency or Laplacian matrix,
#'     depending on \code{model_type}.
#'   \item \strong{coords}: A matrix of \eqn{(x, y)} coordinates.
#'   \item \strong{Q}: The unscaled precision matrix defined by \eqn{W}, \eqn{D}, and \eqn{\gamma}.
#'   \item \strong{phi_true}: The true sampled spatial random effects \eqn{\phi \sim N(0, \tau^2 Q^{-1})}.
#'   \item \strong{eta_true}: The linear predictor \eqn{\eta = X\beta + \phi}.
#'   \item \strong{theta_true}: The relative abundance \eqn{\theta = \exp(\eta)}.
#'   \item \strong{z}: The observed count vector sampled from \eqn{Pois(library\_size * \theta)}.
#' }
#'
#' @note This function requires the \code{Matrix} package for sparse matrix operations.
#'
#' @import Matrix
#' @importFrom Matrix Matrix chol rowSums solve t
#' @importFrom stats arima.sim dist rnorm rpois runif
#'
#' @export
#'
#' @examples
#' set.seed(2026)
#' # Generate data for 1000 cells with a Leroux spatial structure
#' sim_data <- generate_data_one_area(
#'   n_points = 1000,
#'   nb_dist = 0.03,
#'   model_type = "Leroux",
#'   beta_true = c(1, 0, -1),
#'   gamma_true = 0.5,
#'   tau2_true = 1.0,
#'   X_type = "rand_bern"
#' )
generate_data_one_area <- function(n_points,
                                   nb_dist,
                                   model_type,
                                   beta_true,
                                   gamma_true = NULL,
                                   tau2_true = NULL,
                                   X_type = "rand_bern",
                                   ar_gamma = 0.75,
                                   X = NULL,
                                   library_size = NULL) {
  # Model covariance function
  if ("CAR" == model_type) {
    Q_fcn <- Q_matrix_CAR
  } else if ("SAR" == model_type) {
    Q_fcn <- Q_matrix_SAR
  } else if ("Leroux" == model_type) {
    Q_fcn <- Q_matrix_Leroux
  } else {
    stop("Invalid model_type")
  }
  
  # Covariate matrix
  if (is.null(X)) {
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
  }
  
  # Expected value/offset for each point
  if (is.null(library_size)) {
    library_size <- rep(1.0, n_points)
  }
  
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
  
  # D^\{-1\} W for CAR and SAR, D - W for Leroux
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
      coords = coords,
      Q = Q,
      phi_true = phi_true,
      eta_true = eta_true,
      theta_true = theta_true,
      z = z
    )
  )
}


#' Simulate spatial count data via spNNGP
#'
#' Generates synthetic spatial transcriptomics data using a Poisson
#' Generalized Linear Mixed Model (\eqn{GLMM}) with a Sparse Nearest Neighbor
#' Gaussian Process (\eqn{spNNGP}) prior for the random effects. This function
#' approximates a Gaussian Process by using a specified number of nearest
#' neighbors to construct a sparse precision matrix.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param n_points Integer: Number of points/observations in the area.
#' @param nb_dist Numeric: Distance threshold for determining adjacency (if needed
#'   for secondary lattice structures). Coordinates are sampled in \eqn{(0, 1)^2}.
#' @param cov_type Character: The covariance kernel for the Gaussian process.
#'   Options are "Exp" (Exponential), "Mat" (Matern), "Gau" (Gaussian),
#'   and "Sph" (Spherical).
#' @param cov_params Numeric vector: Parameters for the covariance kernel
#'   (typically sigma^2/sill, tau^2/nugget, phi/range, and optionally nu/smoothness).
#' @param nngp_k Integer: The number of nearest neighbors used in the spNNGP
#'   approximation.
#' @param beta_true Numeric vector: True fixed effect coefficients \eqn{\beta}.
#' @param X_type Character: Method for generating the design matrix \eqn{X}.
#'   Options include:
#'   \itemize{
#'     \item \strong{"rand_bern"}: Bernoulli-style covariates \eqn{X_{i, j} \in \{0, 1\}}.
#'     \item \strong{"rand_unif"}: Uniformly distributed \eqn{X_{i, j} \sim U(0, 1)}.
#'     \item \strong{"rand_norm"}: Normally distributed \eqn{X_{i, j} \sim N(0, 1)}.
#'     \item \strong{"ar1"}: \eqn{X} follows an AR(1) process across sorted coordinates.
#'     \item \strong{"ar1_bern"}: AR(1) process transformed into binary values.
#'     \item \strong{"intercept"}: A constant column of 1s.
#'   }
#' @param ar_gamma Numeric: Autocorrelation parameter for \eqn{X} (used if \code{X_type}
#'   is "ar1" or "ar1_bern").
#' @param X Optional matrix: A pre-defined design matrix to use instead of generating one.
#' @param library_size Optional numeric vector: Pre-specified library sizes
#'   (offsets). Defaults to a vector of 1s.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \strong{X}: The generated design matrix of covariates.
#'   \item \strong{W}: The sparse adjacency matrix (if calculated).
#'   \item \strong{D}: The diagonal degree matrix.
#'   \item \strong{eig_list}: Eigenvalues for comparison lattice structures (if calculated).
#'   \item \strong{coords}: A matrix of \eqn{(x, y)} coordinates.
#'   \item \strong{Q}: The sparse precision matrix \eqn{Q = (I - A)^T D^{-1} (I - A)}.
#'   \item \strong{Dinv}: The diagonal variance matrix \eqn{D^{-1}} from the NNGP decomposition.
#'   \item \strong{A}: The sparse lower triangular matrix representing neighbor weights.
#'   \item \strong{phi_true}: The sampled spatial random effects \eqn{\phi \sim N(0, Q^{-1})}.
#'   \item \strong{eta_true}: The linear predictor \eqn{\eta = X\beta + \phi}.
#'   \item \strong{theta_true}: The relative abundance \eqn{\theta = \exp(\eta)}.
#'   \item \strong{z}: The observed count vector sampled from \eqn{Pois(library\_size * \theta)}.
#' }
#'
#' @note This function requires the \code{Matrix} package for handling sparse precision
#'   structures. The spNNGP approach is significantly more memory-efficient than
#'   full Gaussian Processes for large \eqn{n}.
#'
#' @import Matrix
#' @importFrom Matrix Matrix chol Diagonal rowSums solve t
#' @importFrom stats arima.sim dist rnorm rpois runif
#'
#' @export
#'
#' @examples
#' set.seed(2026)
#' # Generate data using a Matern kernel with 20 nearest neighbors
#' sim_data <- generate_data_one_area_spNNGP(
#'   n_points = 1000,
#'   nb_dist = 0.03,
#'   cov_type = "Mat",
#'   cov_params = c(1.0, 0.1, 10, 2), # sill, nugget, range, smoothness
#'   nngp_k = 20,
#'   beta_true = c(1, 0, -1)
#' )
generate_data_one_area_spNNGP <- function(n_points,
                                          nb_dist,
                                          cov_type,
                                          cov_params,
                                          nngp_k,
                                          beta_true,
                                          X_type = "rand_bern",
                                          ar_gamma = 0.75,
                                          X = NULL,
                                          library_size = NULL) {
  # Covariate matrix
  if (is.null(X)) {
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
  }
  
  # Expected value/offset for each point
  if (is.null(library_size)) {
    library_size <- rep(1.0, n_points)
  }
  
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
  if (length(zero_indices) > 0) {
    for (idx in 1:length(zero_indices)) {
      set_idx <- order(d_mat[zero_indices[idx], ])[2]
      W[zero_indices[idx], set_idx] <- 1
      W[set_idx, zero_indices[idx]] <- 1
    }
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
  sp_dist <- sparse_dist_LT(coords, nngp_k)
  close_to_zero_const <- 2.0 * .Machine$double.eps
  
  if ((close_to_zero_const >= cov_params[1]) &&
      (close_to_zero_const >= cov_params[2])) {
    warning("BOTH NUGGET AND SPATIAL VARIANCE ARE ZERO. FAILSAFE: Q IS ZERO AND NOT USED")
    
    Q <- Matrix::Diagonal(n_points, 0)
    A <- Matrix::Diagonal(n_points, 0)
    Dinv <- Matrix::Diagonal(n_points, 0)
  } else {
    # ==========================================
    # --- Precompute Static Neighbor Distance Matrices ---
    n_cols <- ncol(sp_dist)
    nb_dist_list <- vector("list", n_cols)
    
    for (idx in 1:n_cols) {
      keep_idx <- which(!is.na(sp_dist[(1 + nngp_k):(2 * nngp_k), idx]))
      if (length(keep_idx) > 1) {
        nb_dist_list[[idx]] <- as.matrix(stats::dist(coords[sp_dist[nngp_k + keep_idx, idx], , drop = FALSE], diag = TRUE, upper = TRUE))
      } else {
        nb_dist_list[[idx]] <- matrix(0.0, 1, 1)
      }
    }
    # ==========================================
    
    param_est <- nngp_prec_mat(sp_dist, nb_dist_list, cov_type, cov_params)
    Q <- param_est$Q
    A <- param_est$A
    Dinv <- param_est$Dinv
  }
  
  # Spatial random effects
  if ((close_to_zero_const >= cov_params[1]) &&
      (close_to_zero_const >= cov_params[2])) {
    warning("BOTH NUGGET AND SPATIAL VARIANCE ARE ZERO. FAILSAFE: PHI IS ZERO.")
    eta_true <- as.numeric(X %*% beta_true)
    phi_true <- 0 * eta_true
  } else {
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
      coords = coords,
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

#' Simulate spatial counts from a known lattice structure
#'
#' Generates synthetic counts from a Poisson lattice model (CAR, SAR, or Leroux)
#' given a pre-defined adjacency structure and fixed parameters. This is
#' particularly useful for benchmarking against existing spatial graphs or
#' performing power analyses on specific tissue architectures.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param model_type Character: The spatial model for random effects.
#'   Options are "CAR", "SAR", or "Leroux".
#' @param X Matrix: Covariates/design matrix.
#' @param W Matrix: Sparse adjacency matrix representing measurement-to-measurement
#'   connectivity.
#' @param D Matrix: Diagonal degree matrix (row-sums of \code{W}).
#' @param library_size Numeric vector: Scaling factors (offsets) for each
#'   observation.
#' @param tau2_true Numeric: True spatial scale parameter \eqn{\tau^2}.
#' @param gamma_true Numeric: True spatial correlation parameter \eqn{\gamma}.
#' @param beta_true Numeric vector: True fixed effect coefficients \eqn{\beta}.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \strong{z}: The observed count vector sampled from \eqn{Pois(library\_size * \theta)}.
#'   \item \strong{phi_true}: The true sampled spatial random effects
#'     \eqn{\phi \sim N(0, \tau^2 Q^{-1})}.
#'   \item \strong{eta_true}: The linear predictor \eqn{\eta = X\beta + \phi}.
#'   \item \strong{theta_true}: The relative abundance \eqn{\theta = \exp(\eta)}.
#' }
#'
#' @note This function requires the \code{Matrix} package for sparse matrix
#'   operations and Cholesky decomposition.
#'
#' @import Matrix
#' @importFrom Matrix chol solve
#' @importFrom stats rnorm rpois
#'
#' @export
#'
#' @examples
#' set.seed(2026)
#' # Use existing area data to re-sample Poisson counts
#' tau2 <- 1.0
#' gamma <- 0.5
#' beta <- c(1, 0, -1)
#'
#' # Generate base structure
#' base_struct <- generate_data_one_area(1000, 0.03, "Leroux", beta, gamma, tau2, "rand_bern")
#'
#' # Re-sample counts from the same structure
#' sim_counts <- sample_Poisson_lattice(
#'   model_type = "Leroux",
#'   X = base_struct$X,
#'   W = base_struct$W,
#'   D = base_struct$D,
#'   library_size = base_struct$library_size,
#'   tau2_true = tau2,
#'   gamma_true = gamma,
#'   beta_true = beta
#' )
sample_Poisson_lattice <- function(model_type,
                                   X,
                                   W,
                                   D,
                                   library_size,
                                   tau2_true,
                                   gamma_true,
                                   beta_true) {
  # Model covariance function
  if ("CAR" == model_type) {
    Q_fcn <- Q_matrix_CAR
  } else if ("SAR" == model_type) {
    Q_fcn <- Q_matrix_SAR
  } else if ("Leroux" == model_type) {
    Q_fcn <- Q_matrix_Leroux
  } else {
    stop("Invalid model_type")
  }
  
  # Generate inverse precision (unscaled) of random effects
  Q <- Q_fcn(W, D, gamma_true)
  
  # Spatial random effects
  phi_true <- sqrt(tau2_true) * as.numeric(Matrix::solve(Matrix::chol(Q), stats::rnorm(dim(Q)[1])))
  
  # Log Poisson Parameter
  # Add in covariate effect
  eta_true <- as.numeric(X %*% beta_true) + phi_true
  # Get Poisson parameter
  theta_true <- exp(eta_true)
  # Generate data
  z <- stats::rpois(dim(W)[1], theta_true * library_size)
  
  return(list(
    z = z,
    phi_true = phi_true,
    eta_true = eta_true,
    theta_true = theta_true
  ))
}

#' Simulate spatial counts from a known spNNGP structure
#'
#' Generates synthetic counts from a Poisson spNNGP model given a set of spatial
#' coordinates, covariance parameters, and fixed effects. This function is
#' ideal for testing the sensitivity of the spNNGP approximation across
#' different neighborhood sizes \eqn{k} or kernel types without regenerating
#' the underlying spatial layout.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param cov_type Character: The covariance kernel for the Gaussian process.
#'   Options are "Exp" (Exponential), "Mat" (Matern), "Gau" (Gaussian),
#'   and "Sph" (Spherical).
#' @param X Matrix: Covariates/design matrix where rows correspond to observations.
#' @param library_size Numeric vector: Scaling factors (offsets) for each observation.
#' @param coords Matrix: Spatial coordinates \eqn{(x, y)} for each observation.
#'   Must have the same number of rows as \code{X}.
#' @param cov_params Numeric vector: Covariance parameters ordered as
#'   (nugget, sill, range, and optionally smoothness).
#' @param nngp_k Integer: The number of nearest neighbors to use for the
#'   spNNGP precision matrix construction.
#' @param beta_true Numeric vector: True fixed effect coefficients \eqn{\beta}.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \bold{z}: The observed count vector sampled from \eqn{Pois(library\_size * \theta)}.
#'   \item \bold{phi_true}: The true sampled spatial random effects
#'     \eqn{\phi \sim N(0, Q^{-1})}, where \eqn{Q} is the sparse NNGP precision matrix.
#'   \item \bold{eta_true}: The linear predictor \eqn{\eta = X\beta + \phi}.
#'   \item \bold{theta_true}: The relative abundance \eqn{\theta = \exp(\eta)}.
#' }
#'
#' @note This function requires the \code{Matrix} package for sparse Cholesky
#'   decomposition and solving the precision structure.
#'
#' @import Matrix
#' @importFrom Matrix chol solve Diagonal
#' @importFrom stats rnorm rpois dist
#'
#' @export
#'
#' @examples
#' set.seed(2026)
#' # Simulation parameters
#' cov_params <- c(0.1, 1.0, 10, 2) # nugget, sill, range, smoothness
#' beta <- c(1, 0, -1)
#'
#' # Generate initial structure
#' base_data <- generate_data_one_area_spNNGP(
#'   n_points = 1000,
#'   nb_dist = 0.03,
#'   cov_type = "Mat",
#'   cov_params = cov_params,
#'   nngp_k = 20,
#'   beta_true = beta
#' )
#'
#' # Re-sample counts using the known spNNGP kernel
#' sim_counts <- sample_Poisson_spNNGP(
#'   cov_type = "Mat",
#'   X = base_data$X,
#'   library_size = base_data$library_size,
#'   coords = base_data$coords,
#'   cov_params = cov_params,
#'   nngp_k = 20,
#'   beta_true = beta
#' )
sample_Poisson_spNNGP <- function(cov_type,
                                  X,
                                  library_size,
                                  coords,
                                  cov_params,
                                  nngp_k,
                                  beta_true) {
  # Generate inverse precision (unscaled) of random effects
  sp_dist <- sparse_dist_LT(coords, nngp_k)
  close_to_zero_const <- 2.0 * .Machine$double.eps
  
  if ((close_to_zero_const >= cov_params[1]) &&
      (close_to_zero_const >= cov_params[2])) {
    warning("BOTH NUGGET AND SPATIAL VARIANCE ARE ZERO. FAILSAFE: Q IS ZERO AND NOT USED")
    Q <- Matrix::Diagonal(nrow(coords), 0)
  } else {
    # ==========================================
    # --- Precompute Static Neighbor Distance Matrices ---
    n_cols <- ncol(sp_dist)
    nb_dist <- vector("list", n_cols)
    
    for (idx in 1:n_cols) {
      keep_idx <- which(!is.na(sp_dist[(1 + nngp_k):(2 * nngp_k), idx]))
      if (length(keep_idx) > 1) {
        nb_dist[[idx]] <- as.matrix(stats::dist(coords[sp_dist[nngp_k + keep_idx, idx], , drop = FALSE], diag = TRUE, upper = TRUE))
      } else {
        nb_dist[[idx]] <- matrix(0.0, 1, 1)
      }
    }
    # ==========================================
    
    param_est <- nngp_prec_mat(sp_dist, nb_dist, cov_type, cov_params)
    Q <- param_est$Q
  }
  
  # Spatial random effects
  if ((close_to_zero_const >= cov_params[1]) &&
      (close_to_zero_const >= cov_params[2])) {
    warning("BOTH NUGGET AND SPATIAL VARIANCE ARE ZERO. FAILSAFE: PHI IS ZERO.")
    eta_true <- as.numeric(X %*% beta_true)
    phi_true <- 0 * eta_true
  } else {
    phi_true <- as.numeric(Matrix::solve(Matrix::chol(Q), stats::rnorm(nrow(coords))))
    # Add in covariate effect
    eta_true <- as.numeric(X %*% beta_true) + phi_true
  }
  
  # Get Poisson parameter
  theta_true <- exp(eta_true)
  # Generate data
  z <- stats::rpois(nrow(coords), theta_true * library_size)
  
  return(list(
    z = z,
    phi_true = phi_true,
    eta_true = eta_true,
    theta_true = theta_true
  ))
}

#' Generate synthetic multi-sample spatial count data
#'
#' Replaces the observed counts in a \code{TESSERAData_obj} with synthetic counts
#' generated from specified true parameters. This utility is designed for
#' benchmarking and power analyses, allowing users to simulate data over
#' existing tissue architectures and experimental designs.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param TESSERAData_obj Object containing the experimental structure,
#'   typically created by the \code{prepData} method.
#' @param gene_list Character vector: The names of the genes/measurements to
#'   simulate. These must correspond to row names in the original count data.
#' @param data_gen_model Character: The spatial model to use for simulations.
#'   Options include "CAR", "SAR", "Leroux" (Lattice models), or "spNNGP".
#' @param tau2_true Spatial scale parameters \eqn{\tau^2}. Provide a numeric vector
#'   (for a single gene) or a matrix of dimensions (genes x samples) (for
#'   multiple genes). Required only for Lattice models.
#'   \strong{Note:} \eqn{\tau^2 \ge 0} is required.
#' @param gamma_true Spatial correlation parameters \eqn{\gamma}. Provide a numeric
#'   vector (single gene) or a matrix of dimensions (genes x samples).
#'   Required only for Lattice models.
#'   \itemize{
#'     \item \strong{CAR/SAR}: \eqn{-1 < \gamma < 1} required.
#'     \item \strong{Leroux}: \eqn{0 \le \gamma < 1} required.
#'   }
#' @param cov_params Spatial covariance parameters for spNNGP. Provide a matrix
#'   (samples x 3/4 parameters) for a single gene, or an array
#'   (genes x samples x 3/4 parameters) for multiple genes.
#'   Parameters must be ordered as (nugget, sill, range, and optional smoothness).
#' @param cov_type Character: The spatial correlation kernel for spNNGP.
#'   Options are "Exp", "Mat", "Gau", and "Sph".
#' @param nngp_k Integer: The number of nearest neighbors to use for the
#'   spNNGP kernel approximation.
#' @param beta_true True fixed effect coefficients \eqn{\beta}. Provide a numeric
#'   vector (single gene) or a matrix of dimensions (genes x covariates).
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \strong{new_TESSERAData_obj}: A \code{TESSERAData_obj} where the
#'     \code{counts_list} has been replaced with synthetic values.
#'   \item \strong{synthetic_count_summary}: A data frame summarizing the
#'     generated counts and computing similarity metrics between the
#'     synthetic and original data.
#' }
#'
#' @importFrom dplyr bind_rows
#'
#' @export
#'
#' @examples
#' # Locate the prepped TESSERAData object in inst/extdata
#' rds_path <- system.file("extdata", "example_prepData.rds", package = "TESSERA")
#' # Load the TESSERAData object
#' TESSERA_data <- readRDS(rds_path)
#'
#' #' # Locate the saved model results in inst/extdata
#' rds_path <- system.file("extdata", "example_TESSERA_out_Leroux.rds", package = "TESSERA")
#' # Load the results object
#' TESSERA_out_Leroux <- readRDS(rds_path)
#'
#' # Run the synthetic data generation
#' TESSERA_resampled_data <- prep_synth_data(
#'   TESSERAData_obj = TESSERA_data,
#'   gene_list = "example",
#'   data_gen_model = "Leroux",
#'   tau2_true = TESSERA_out_Leroux$tau2_hat,
#'   gamma_true = TESSERA_out_Leroux$gamma_hat,
#'   beta_true = TESSERA_out_Leroux$beta_hat
#' )$new_TESSERAData_obj
prep_synth_data <- function(TESSERAData_obj,
                            gene_list,
                            data_gen_model,
                            tau2_true = NULL,
                            gamma_true = NULL,
                            cov_params = NULL,
                            cov_type = NULL,
                            nngp_k = NULL,
                            beta_true = NULL) {
  # Make sure object is correct/usable
  if ("spNNGP" == data_gen_model) {
    check_inputs_TESSERA_spNNGP(TESSERAData_obj)
  } else if (("CAR" == data_gen_model) ||
             ("SAR" == data_gen_model) ||
             ("Leroux" == data_gen_model)) {
    check_inputs_TESSERA(TESSERAData_obj)
  } else {
    stop("Invalid data_gen_model.")
  }
  
  
  # Handle case of only one gene
  if (!is.vector(gene_list) && !is.list(gene_list)) {
    gene_list <- c(gene_list)
  }
  if (is.list(gene_list)) {
    gene_list <- c(unlist(gene_list))
  }
  
  # Convert inputs to one-row matrices (Lattice models) or arrays
  # Also check dimensions
  
  # Check beta
  if (!is.matrix(beta_true)) {
    beta_true <- matrix(data = beta_true,
                        nrow = 1,
                        ncol = length(beta_true))
  }
  stopifnot(length(gene_list) == nrow(beta_true))
  stopifnot(ncol(TESSERAData_obj$X_list[[1]]) == ncol(beta_true))
  
  # Check gamma, tau^2 for Lattice models
  if (("CAR" == data_gen_model) ||
      ("SAR" == data_gen_model) || ("Leroux" == data_gen_model)) {
    if (!is.matrix(tau2_true)) {
      tau2_true <- matrix(data = tau2_true,
                          nrow = 1,
                          ncol = length(tau2_true))
    }
    if (!is.matrix(gamma_true)) {
      gamma_true <- matrix(data = gamma_true,
                           nrow = 1,
                           ncol = length(gamma_true))
    }
    
    # Check dimensions of inputs
    stopifnot(length(gene_list) == nrow(tau2_true))
    stopifnot(length(gene_list) == nrow(gamma_true))
    stopifnot(length(TESSERAData_obj$counts_list) == ncol(tau2_true))
    stopifnot(length(TESSERAData_obj$counts_list) == ncol(gamma_true))
  } else if ("spNNGP" == data_gen_model) {
    # genes X samples X parameters
    if ((1 == length(dim(cov_params))) ||
        (0 == length(dim(cov_params)))) {
      # Only one gene, one sample passed in
      tmp <- array(data = 0,
                   dim = c(
                     length(gene_list),
                     length(TESSERAData_obj$X_list),
                     length(cov_params)
                   ))
      tmp[1, 1, ] <- cov_params
      cov_params <- tmp
    } else if (2 == length(dim(cov_params))) {
      # Only one gene, but multiple samples passed in
      tmp <- array(data = 0,
                   dim = c(
                     length(gene_list),
                     length(TESSERAData_obj$X_list),
                     ncol(cov_params)
                   ))
      tmp[1, , ] <- cov_params
      cov_params <- tmp
    }
    # Check dimensions of inputs
    stopifnot(length(gene_list) == dim(cov_params)[1])
    stopifnot(length(TESSERAData_obj$counts_list) == dim(cov_params)[2])
    stopifnot(3 == length(dim(cov_params)))
  }
  
  # Synthetic counts list
  synth_counts_list <- list()
  # Initialize synthetic counts list with matrices for each sample
  for (s_idx in 1:length(TESSERAData_obj$counts_list)) {
    synth_counts_list[[s_idx]] <- matrix(0,
                                         nrow = length(gene_list),
                                         ncol = ncol(TESSERAData_obj$counts_list[[s_idx]]))
    rownames(synth_counts_list[[s_idx]]) <- gene_list
    colnames(synth_counts_list[[s_idx]]) <- colnames(TESSERAData_obj$counts_list[[s_idx]])
  }
  names(synth_counts_list) <- names(TESSERAData_obj$counts_list)
  
  # Sample from a Poisson lattice model to obtain new count matrices for each sample
  summary_df <- list()
  for (g_idx in 1:length(gene_list)) {
    gene <- gene_list[g_idx]
    
    for (s_idx in 1:length(TESSERAData_obj$counts_list)) {
      if (("CAR" == data_gen_model) ||
          ("SAR" == data_gen_model) ||
          ("Leroux" == data_gen_model)) {
        synth_counts_list[[s_idx]][g_idx, ] <- sample_Poisson_lattice(
          model_type = data_gen_model,
          X = TESSERAData_obj$X_list[[s_idx]],
          W = TESSERAData_obj$W_list[[s_idx]],
          D = TESSERAData_obj$D_list[[s_idx]],
          library_size = TESSERAData_obj$library_size_list[[s_idx]],
          tau2_true = tau2_true[g_idx, s_idx],
          gamma_true = gamma_true[g_idx, s_idx],
          beta_true = beta_true[g_idx, ]
        )$z
      } else if ("spNNGP" == data_gen_model) {
        synth_counts_list[[s_idx]][g_idx, ] <- sample_Poisson_spNNGP(
          cov_type = cov_type,
          X = TESSERAData_obj$X_list[[s_idx]],
          library_size = TESSERAData_obj$library_size_list[[s_idx]],
          coords = TESSERAData_obj$coords_list[[s_idx]],
          cov_params = cov_params[g_idx, s_idx, ],
          nngp_k = nngp_k,
          beta_true = beta_true[g_idx, ]
        )$z
      }
    }
    
    # Store summary about the generated counts/comparison with real data
    synth_summary <- sapply(1:length(TESSERAData_obj$counts_list), function (x) {
      x_tmp <- as.vector(synth_counts_list[[x]][g_idx, ])
      y_tmp <- as.vector(TESSERAData_obj$counts_list[[x]][g_idx, ])
      
      # Create PMFs for both real/synthetic data
      tX <- as.data.frame(table(x_tmp))
      colnames(tX) <- c("Val", "X")
      tY <- as.data.frame(table(y_tmp))
      colnames(tY) <- c("Val", "Y")
      tXY <- merge(tX, tY, all = TRUE)
      tXY[is.na(tXY)] <- 0
      tXY$X <- tXY$X / sum(tXY$X)
      tXY$Y <- tXY$Y / sum(tXY$Y)
      
      # Summary information
      s_X <- summary(as.vector(synth_counts_list[[x]][g_idx, ]))
      names(s_X) <- paste0("synth_", names(s_X))
      s_Y <- summary(as.vector(TESSERAData_obj$counts_list[[x]][g_idx, ]))
      names(s_Y) <- paste0("orig_", names(s_Y))
      return(Reduce(cbind, list(
        data.frame(
          gene = gene,
          sample = names(TESSERAData_obj$counts_list)[x],
          "TV" = sum(abs(tXY$X - tXY$Y)) / 2,
          # Total variation distance [0, 1]
          "KS" = max(abs(tXY$X - tXY$Y)),
          # KS statistic [0, 1]
          "COR" = cor(x_tmp, y_tmp, method = "spearman"),
          # Correlation [-1, 1]
          "NormMSE" = mean((x_tmp - y_tmp)^2) / mean(y_tmp^2)
          # Normalized MSE [0, \infty),
        ),
        data.frame(as.list(s_X)),
        data.frame(as.list(s_Y))
      )))
    })
    synth_summary <- t(synth_summary)
    summary_df[[1 + length(summary_df)]] <- data.frame(synth_summary)
  }
  if (1 < length(gene_list)) {
    summary_df <- dplyr::bind_rows(summary_df)
  }
  else {
    summary_df <- summary_df[[1]]
  }
  
  # Create new TESSERAData object
  new_TESSERAData_obj <- TESSERAData_obj
  new_TESSERAData_obj$counts_list <- synth_counts_list
  
  return(
    list(
      new_TESSERAData_obj = new_TESSERAData_obj,
      synthetic_count_summary = summary_df
    )
  )
}
