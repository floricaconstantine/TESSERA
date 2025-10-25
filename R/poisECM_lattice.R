## Wrapper functions to fit model.
# Dependencies in file: Matrix.
# Dependencies: Functions from utils.R, models.R, E_step.R, M_step.R.
# Dependences from functions in other files not listed: pracma, sp, gstat.
# Rcpp dependencies: calc_moran.cpp.


#' Fit multi-area Poisson spatial generalized linear model.
#'  Allows for CAR, SAR, or Leroux random effects.
#'  Fits a common set of fixed effects across areas.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param poisECMData_obj Object containing data.
#'  Created by the prepData method.
#' @param gene_name Which gene/measurement to fit.
#'  I.e., Which row in the count data matrix to fit.
#' @param model_type Which model to fit for the random effects.
#'    "CAR", "SAR", and "Leroux" are the valid options.
#' @param em_iters Number of E(C)M iterations to run.
#' @param opt_iters Number of inner M-step iterations to run.
#'    The M-step is technically a CM step, i.e., a conditional maximization.
#'    That is, we perform block coordinate ascent wherein we maximize the
#'    expected likelihood in tau^2 and then in gamma for each area, and then
#'    in beta, where each maximization holds all other parameters constant.
#'    This approach leads to a slower convergence in terms of the number of
#'    iterations (em_iters) but faster computation in each iteration.
#' @param em_min_iters Minimum number of E(C)M iterations to run.
#' @param em_tol Tolerance for early stopping the E(C)M algorithm.
#' @param em_stopping How to decide to stop early.
#'    NULL: Don't stop early.
#'    "abs_loglike": Difference in absolute data log likelihood, |old - new|.
#'    "rel_loglike": Relative difference in data log likelihood, |old - new|/|old|.
#'      Data log likelihood taken across alll areas.
#'    "abs_beta_norm": Absolute L2 norm of difference in beta, ||old - new||_2.
#'    "rel_beta_norm": Relative L2 norm of difference in beta,
#'      ||old - new||_2 / ||old||_2.
#' @param beta_init How to initialize beta.
#'    "glm": Stack up all the observed counts z and covariate matrices X and
#'      fit a Poisson GLM z ~ 0 + X using base R's glm function.
#'    "random": Initialize with standard normal random variables.
#'    A vector or matrix: Initialize with provided value.
#' @param gamma_init How to initialize beta.
#'    "moran": Compute absolute Moran's I of log(z + 1/2) - X beta
#'    "random": Initialize with standard uniform random variables.
#'    A number, vector, or matrix: Initialize with provided values.
#' @param tau2_init How to initialize beta.
#'    "lognormal": Approximation based on the log-normality of the Poisson parameter.
#'    "var": Compute variance of log(z + 1/2) - X beta
#'    "random": Initialize with standard normal random variables.
#'    A number, vector, or matrix: Initialize with provided values.
#' @param verbose Boolean of whether to print out updates.
#' @param dense_matrices Boolean of whether to treat Q as a dense matrix when
#'    computing certain quantities in the E-step. This will increase memory usage
#'    (e.g., ~2 -> 8 GB for 3k cells per sample and 20 neighbors), but will slightly
#'    speed up computation (by roughly 10-20%).
#'
#' @return A list comprised of the following:
#' @return beta_hat: Estimated coefficients.
#' @return gamma_hat: Estimated spatial correlation parameters.
#' @return tau2_hat: Estimated spatial scale parameters.
#' @returns predictions: Predicted values.
#'    Vector of all values (corresponds to z_list stacked together).
#' @returns eta_hat: Estimated random effects.
#'    Vector of all values (corresponds to z_list stacked together).
#' @returns theta_hat: Estimated Poisson parameters.
#'    Vector of all values (corresponds to z_list stacked together).
#' @returns phi_hat: Estimated spatial random effects.
#'    Vector of all values (corresponds to z_list stacked together).
#' @returns gamma_tracker: Estimated correlation parameters.
#'    areas x EM iterations matrix---history across iterations.
#' @returns tau2_tracker: Estimated scale parameters.
#'    areas x EM iterations matrix---history across iterations.
#' @returns beta_tracker: Estimated fixed effects.
#'    coefficients x EM iterations matrix---history across iterations.
#' @returns R2_tracker: Squared correlation between predictions and observations.
#'    areas x EM Iterations matrix---history across iterations.
#' @returns MSE_tracker: Mean-Square Error between predictions and observations.
#'    areas x EM Iterations matrix---history across iterations.
#' @returns data_log_like_tracker: Log likelihood of observations.
#'    areas x EM Iterations matrix---history across iterations.
#' @returns expected_log_like_tracker: Expected log likelihood at current parameters.
#'    areas x EM Iterations matrix---history across iterations.
#' @returns resid_moran_nb: Moran's I for each area computed using neighbor
#'   adjacency as weights.
#'   areas x EM Iterations matrix---history across iterations.
#' @returns resid_moran: Moran's I for each area computed using coordinates.
#'   areas x EM Iterations x (value, expectation, sd) array.
#'   NULL unless coordinates are supplied.
#' @returns start_idx_list: Indices where each area's values start in predictions, etc.
#'    E.g., if area 1 has 100 points and area 2 has 50, we would have
#'    c(1, 101, 151, ...).
#' @returns eig_val_list: List of eigenvalues, depending on model type.
#'    ONLY IF NOT PASSED IN.
#' @returns time: Total time taken by function.
#' @returns beta_neghessian: Negative Hessian of final value of beta.
#'    Matrix.
#' @returns tau2_neghessian: Negative Hessian of final value of tau^2.
#'    Vector of values.
#' @returns gamma_neghessian: Negative Hessian of final value of gamma.
#'    Vector of values.
#' @returns run_settings: A list with the parameter settings used to run the algorithm.
#' @returns performanceSummary: A dataframe with summary statistics for each sample.
#'  See [summarizePoisECMPerformance()] for more details.
#'
#' @references Meng, Xiao-Li, and Donald B. Rubin.
#'                "Maximum likelihood estimation via the ECM algorithm: A general framework."
#'                Biometrika 80.2 (1993): 267-278.
#'
#' @note D must be invertible, as must W, so make sure that every point has
#'      at least one neighbor, i.e., that there are no isolated points,
#'      for at least the CAR and SAR models.
#' @note Requires the Matrix library.
#'
#' @import Matrix
#' @importFrom stats coef
#' @importFrom stats cor
#' @importFrom stats dpois
#' @importFrom stats poisson
#' @importFrom stats predict
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats var
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @useDynLib poisECM
#'
#' @export
poisECM_lattice <- function(poisECMData_obj,
                            gene_name,
                            model_type = "CAR",
                            em_iters = 200,
                            opt_iters = 10,
                            em_min_iters = 10,
                            em_tol = 1e-3,
                            em_stopping = NULL,
                            beta_init = "glm",
                            gamma_init = "moran",
                            tau2_init = "var",
                            verbose = TRUE,
                            dense_matrices = FALSE) {
  # Start clock
  t0_EM = Sys.time()

  # Check inputs
  checkInputsPoisECM(poisECMData_obj)

  # Extract gene of interest and associated counts
  gene_idx <- which(rownames(poisECMData_obj$counts_list[[1]]) == gene_name)
  z_list <- lapply(poisECMData_obj$counts_list, function (x) {
    x[gene_idx, ]
  })

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

  # gamma optimization function
  if ("CAR" == model_type) {
    gamma_M_fcn <- M_step_gamma_CAR
  } else if ("SAR" == model_type) {
    gamma_M_fcn = M_step_gamma_SAR
  } else if ("Leroux" == model_type) {
    gamma_M_fcn <- M_step_gamma_Leroux
  } else {
    stop("Invalid model_type")
  }

  # Number of areas
  n_areas <- length(poisECMData_obj$W_list)
  # Total number of points
  n_total_points <- 0
  # Starting points in larger matrices
  start_idx_list <- rep(NA, length(z_list))
  for (idx in 1:length(poisECMData_obj$X_list)) {
    n_total_points <- n_total_points + nrow(poisECMData_obj$X_list[[idx]])
    if (1 == idx) {
      start_idx_list[idx] <- 1
    } else {
      start_idx_list[idx] <- start_idx_list[idx - 1] + length(z_list[[idx - 1]])
    }
  }
  # Dimension of coefficient vector
  beta_dim = ncol(poisECMData_obj$X_list[[1]])

  # Compute eigenvalues if needed
  if (is.null(poisECMData_obj$eig_CS_list) &&
      is.null(poisECMData_obj$eig_L_list)) {
    cat("Eigenvalues not supplied: Computing now.")
    eig_val_list <- list()
    for (area_idx in 1:n_areas) {
      # D^{-1} W for CAR and SAR, D - W for Leroux
      if (("CAR" == model_type) || ("SAR" == model_type)) {
        eig_val_list[[area_idx]] <- Re(eigen(
          Matrix::solve(poisECMData_obj$D_list[[area_idx]], poisECMData_obj$W_list[[area_idx]]),
          FALSE,
          only.values = TRUE
        )$values)
      } else if ("Leroux" == model_type) {
        eig_val_list[[area_idx]] <- Re(
          eigen(
            poisECMData_obj$D_list[[area_idx]] - poisECMData_obj$W_list[[area_idx]],
            TRUE,
            only.values = TRUE
          )$values
        )
      }
    }
  } else {
    # Eigenvalues
    if ((model_type == "CAR") || (model_type == "SAR")) {
      eig_val_list <- poisECMData_obj$eig_CS_list
    } else if (model_type == "Leroux") {
      eig_val_list <- poisECMData_obj$eig_L_list
    }
  }

  # Store parameter estimates and track
  gamma_tracker <- matrix(data = NA,
                          nrow = n_areas,
                          ncol = em_iters + 1)
  tau2_tracker <- matrix(data = NA,
                         nrow = n_areas,
                         ncol = em_iters + 1)
  beta_tracker <- matrix(data = NA,
                         nrow = beta_dim,
                         ncol = em_iters + 1)
  fit_tracker <- matrix(data = NA,
                        nrow = n_total_points,
                        ncol = em_iters)
  eta_tracker <- matrix(data = NA,
                        nrow = n_total_points,
                        ncol = em_iters)
  theta_tracker <- matrix(data = NA,
                          nrow = n_total_points,
                          ncol = em_iters)

  # Track performance
  R2_tracker <- matrix(data = NA,
                       nrow = n_areas,
                       ncol = em_iters)
  MSE_tracker <- matrix(data = NA,
                        nrow = n_areas,
                        ncol = em_iters)
  data_log_like_tracker <- matrix(data = NA,
                                  nrow = n_areas,
                                  ncol = em_iters)
  expected_log_like_tracker <- matrix(data = NA,
                                      nrow = n_areas,
                                      ncol = em_iters)
  resid_moran <- array(data = NA, dim = c(n_areas, em_iters, 3))
  resid_moran_nb <- matrix(data = NA,
                           nrow = n_areas,
                           ncol = em_iters)

  # Initialize parameters: beta
  if (is.character(beta_init) && ("random" == beta_init)) {
    beta_tracker[, 1] <- stats::rnorm(beta_dim)
    cat("Random initialization for beta.", "\n")
  } else if (is.character(beta_init) && ("glm" == beta_init)) {
    # Let's be intelligent: initialize with a basic GLM

    # Stack into a single vector/matrix
    z_vec <- Reduce(c, z_list)
    lib_vec <- Reduce(c, poisECMData_obj$library_size_list)
    X_mat <- Reduce(rbind, poisECMData_obj$X_list)
    # Fit GLM
    beta_tmp <- as.vector(stats::coef(
      # stats::glm(z_vec / lib_vec ~ 0 + X_mat, family = stats::poisson())
      stats::glm(
        z_vec ~ 0 + X_mat,
        family = stats::poisson(),
        offset = log(lib_vec)
      )
    ))
    beta_tmp[is.nan(beta_tmp)] <- 0
    beta_tmp[is.infinite(beta_tmp)] <- 0
    beta_tracker[, 1] <- beta_tmp

    # Memory
    rm(z_vec)
    rm(X_mat)

    cat("GLM initialization for beta.", "\n")
  } else if (is.numeric(beta_init) &&
             (is.vector(beta_init) || is.matrix(beta_init))) {
    # Pass in a value
    beta_tracker[, 1] <- as.vector(beta_init)

    cat("Pre-defined initialization for beta.", "\n")
  } else {
    stop("Invalid initialization for beta.")
  }
  # Handle NA
  beta_tracker[is.nan(beta_tracker[, 1]), 1] <- 0.0
  beta_tracker[is.na(beta_tracker[, 1]), 1] <- 0.0
  cat("Initial beta", beta_tracker[, 1], "\n")

  # Initialize parameters: gamma
  if (is.character(gamma_init) && ("moran" == gamma_init)) {
    for (area_idx in 1:n_areas) {
      tmp_z <- log((0.5 + z_list[[area_idx]]) / poisECMData_obj$library_size_list[[area_idx]]) - poisECMData_obj$X_list[[area_idx]] %*% beta_tracker[, 1]
      gamma_tracker[area_idx, 1] <- abs(moran_I_nb(tmp_z, poisECMData_obj$W_list[[area_idx]]))
    }
    cat("Moran initialization for gamma", "\n")
  } else if (is.character(gamma_init) && ("random" == gamma_init)) {
    gamma_tracker[, 1] <- stats::runif(n_areas)
    cat("Random initialization for gamma", "\n")
  } else if (is.numeric(gamma_init)) {
    gamma_tracker[, 1] <- gamma_init
    cat("Predefined initialization for gamma", "\n")
  } else {
    stop("Invalid initialization for gamma.")
  }
  # Handle zero counts/NaN values
  gamma_tracker[is.nan(gamma_tracker[, 1]), 1] <- 0.0
  gamma_tracker[is.na(gamma_tracker[, 1]), 1] <- 0.0

  cat("Initial gamma", gamma_tracker[, 1], "\n")



  # Initialize parameters: tau^2
  if (is.character(tau2_init) && ("lognormal" == tau2_init)) {
    for (area_idx in 1:n_areas) {
      glm_out <- stats::glm((z_list[[area_idx]] / poisECMData_obj$library_size_list[[area_idx]]) ~ 0 + poisECMData_obj$X_list[[area_idx]],
                            family = stats::poisson()
      )
      tau2_tracker[area_idx, 1] <- max(
        1e-6,
        2 * log(mean(
          (z_list[[area_idx]] / poisECMData_obj$library_size_list[[area_idx]])
          / stats::predict(glm_out, type = "response")
        )) #- 2
      )
    }
    cat("Log-Normal initialization for tau^2", "\n")
  } else if (is.character(tau2_init) && ("var" == tau2_init)) {
    for (area_idx in 1:n_areas) {
      tmp_z <- log((0.5 + z_list[[area_idx]]) / poisECMData_obj$library_size_list[[area_idx]]) - poisECMData_obj$X_list[[area_idx]] %*% beta_tracker[, 1]
      tau2_tracker[area_idx, 1] <- stats::var(tmp_z)
    }
    cat("Variance initialization for tau^2", "\n")
  } else if (is.character(tau2_init) && ("random" == tau2_init)) {
    tau2_tracker[, 1] <- abs(stats::rnorm(n_areas))
    cat("Random initialization for tau^2", "\n")
  } else if (is.numeric(tau2_init)) {
    tau2_tracker[, 1] <- tau2_init
    cat("Predefined initialization for tau^2", "\n")
  } else {
    stop("Invalid initialization for tau^2.")
  }
  # Handle zero counts/NaN values
  tau2_tracker[is.nan(tau2_tracker[, 1]), 1] <- 0.0
  tau2_tracker[is.na(tau2_tracker[, 1]), 1] <- 0.0
  tau2_tracker[tau2_tracker[, 1] == 0.0, 1] <- 1e-3
  cat("Initial tau2", tau2_tracker[, 1], "\n")

  # Also initialize dependency Q as a function of W, D, and gamma
  Q_hat_list <- list()
  for (area_idx in 1:n_areas) {
    Q_hat_list[[area_idx]] <- Q_fcn(poisECMData_obj$W_list[[area_idx]],
                                    poisECMData_obj$D_list[[area_idx]],
                                    gamma_tracker[area_idx, 1])
  }

  # Loop over EM iterations
  for (em_idx in 1:em_iters) {
    if (verbose || (0 == (em_idx %% 100))) {
      cat(
        paste(
          "Start of EM Iteration",
          em_idx,
          "of",
          em_iters,
          ";",
          Sys.time() - t0_EM,
          "Elapsed"
        ),
        "\n"
      )
    }



    # Run E-Step: Get covariance and mean of eta
    Vhat_list <- list()
    eta_hat_list <- list()
    for (area_idx in 1:n_areas) {
      if (dense_matrices) {
        Vhat_list[[area_idx]] <- Vhat_list[[area_idx]] <- E_step_Vhat(as.matrix(Q_hat_list[[area_idx]]), tau2_tracker[area_idx, em_idx], z_list[[area_idx]])
      } else {
        Vhat_list[[area_idx]] <- E_step_Vhat(Q_hat_list[[area_idx]], tau2_tracker[area_idx, em_idx], z_list[[area_idx]])
      }

      eta_hat_list[[area_idx]] <- E_step_etahat(
        Vhat_list[[area_idx]],
        Q_hat_list[[area_idx]],
        tau2_tracker[area_idx, em_idx],
        beta_tracker[, em_idx],
        poisECMData_obj$X_list[[area_idx]],
        z_list[[area_idx]],
        poisECMData_obj$library_size_list[[area_idx]]
      )
      # Get a few more things out of the E-step: theta and predictions
      theta_hat <- as.numeric(E_step_thetahat(Vhat_list[[area_idx]], eta_hat_list[[area_idx]]))
      z_hat <- as.numeric(E_step_predict(theta_hat, poisECMData_obj$library_size_list[[area_idx]]))

      # Store stuff
      fit_tracker[start_idx_list[area_idx]:(start_idx_list[area_idx] + length(z_hat) - 1), em_idx] <- z_hat
      eta_tracker[start_idx_list[area_idx]:(start_idx_list[area_idx] + length(eta_hat_list[[area_idx]]) - 1), em_idx] <- as.numeric(eta_hat_list[[area_idx]])
      theta_tracker[start_idx_list[area_idx]:(start_idx_list[area_idx] + length(theta_hat) - 1), em_idx] <- theta_hat

      # Performance
      R2_tracker[area_idx, em_idx] <- stats::cor(z_hat, z_list[[area_idx]])^2
      MSE_tracker[area_idx, em_idx] <- mean(abs(z_hat - z_list[[area_idx]])^2)
      # Observed, expected, sd
      if (!(
        is.null(poisECMData_obj$coords_list[[area_idx]][, 1]) ||
        is.null(poisECMData_obj$coords_list[[area_idx]][, 2])
      )) {
        resid_moran[area_idx, em_idx, ] <- calc_moran(
          z_hat - z_list[[area_idx]],
          poisECMData_obj$coords_list[[area_idx]][, 1],
          poisECMData_obj$coords_list[[area_idx]][, 2]
        )
      }
      resid_moran_nb[area_idx, em_idx] <- moran_I_nb(z_hat - z_list[[area_idx]], poisECMData_obj$W_list[[area_idx]])

      # Data log likelihood
      # data_log_like_tracker[area_idx, em_idx] <- poisson_loglike(z_list[[area_idx]], theta_hat * library_size_list[[area_idx]])
      data_log_like_tracker[area_idx, em_idx] <- sum(
        stats::dpois(
          round(z_list[[area_idx]]),
          theta_hat * poisECMData_obj$library_size_list[[area_idx]],
          log = TRUE
        ),
        na.rm = TRUE
      )

      # Expected log likelihood
      expected_log_like_tracker[area_idx, em_idx] <- expected_loglike(
        Vhat_list[[area_idx]],
        eta_hat_list[[area_idx]],
        Q_hat_list[[area_idx]],
        gamma_tracker[area_idx, em_idx],
        tau2_tracker[area_idx, em_idx],
        beta_tracker[, em_idx],
        poisECMData_obj$X_list[[area_idx]],
        poisECMData_obj$W_list[[area_idx]],
        poisECMData_obj$D_list[[area_idx]],
        eig_val_list[[area_idx]],
        model_type
      )
    }

    # Run (C)M-Step: Optimize
    for (opt_idx in 1:opt_iters) {
      # Initialize with current values
      if (1 == opt_idx) {
        gamma_tracker[, 1 + em_idx] <- gamma_tracker[, em_idx]
        tau2_tracker[, 1 + em_idx] <- tau2_tracker[, em_idx]
        beta_tracker[, 1 + em_idx] <- beta_tracker[, em_idx]
      }

      for (area_idx in 1:n_areas) {
        # Optimize in tau^2
        tau2_tracker[area_idx, 1 + em_idx] <- M_step_tau2(
          Vhat_list[[area_idx]],
          eta_hat_list[[area_idx]],
          Q_hat_list[[area_idx]],
          beta_tracker[, 1 + em_idx],
          poisECMData_obj$X_list[[area_idx]]
        )

        # Optimize in gamma
        gamma_out <- gamma_M_fcn(
          Vhat_list[[area_idx]],
          eta_hat_list[[area_idx]],
          tau2_tracker[area_idx, 1 + em_idx],
          beta_tracker[, 1 + em_idx],
          poisECMData_obj$X_list[[area_idx]],
          poisECMData_obj$W_list[[area_idx]],
          poisECMData_obj$D_list[[area_idx]],
          eig_val_list[[area_idx]],
          gamma_tracker[area_idx, 1 + em_idx]
        )
        # Handle optimization going off the rails
        if (("CAR" == model_type) | ("SAR" == model_type)) {
          if (-1 > gamma_out$gamma_hat) {
            warning(
              paste0(
                "Iteration ",
                em_idx,
                ", Inner iteration ",
                opt_idx,
                ", Area ",
                area_idx,
                " gamma is less than -1, correcting"
              )
            )
            gamma_out$gamma_hat <- max(-1, gamma_out$gamma_hat)
          }
          if (1 < gamma_out$gamma_hat) {
            warning(
              paste0(
                "Iteration ",
                em_idx,
                ", Inner iteration ",
                opt_idx,
                ", Area ",
                area_idx,
                " gamma is larger than 1, correcting"
              )
            )
            gamma_out$gamma_hat <- min(1, gamma_out$gamma_hat)
          }
        } else if ("Leroux" == model_type) {
          if (0 > gamma_out$gamma_hat) {
            warning(
              paste0(
                "Iteration ",
                em_idx,
                ", Inner iteration ",
                opt_idx,
                ", Area ",
                area_idx,
                " gamma is less than 0, correcting"
              )
            )
            gamma_out$gamma_hat <- max(0, gamma_out$gamma_hat)
          }
          if (1 < gamma_out$gamma_hat) {
            warning(
              paste0(
                "Iteration ",
                em_idx,
                ", Inner iteration ",
                opt_idx,
                ", Area ",
                area_idx,
                " gamma is larger than 1, correcting"
              )
            )
            gamma_out$gamma_hat <- min(1, gamma_out$gamma_hat)
          }
        }
        gamma_tracker[area_idx, 1 + em_idx] <- gamma_out$gamma_hat

        # Update precision matrix Q
        Q_hat_list[[area_idx]] <- Q_fcn(poisECMData_obj$W_list[[area_idx]],
                                        poisECMData_obj$D_list[[area_idx]],
                                        gamma_tracker[area_idx, em_idx + 1])
      }

      # Optimize in beta
      beta_tracker[, 1 + em_idx] <- M_step_beta(eta_hat_list,
                                                Q_hat_list,
                                                tau2_tracker[, 1 + em_idx],
                                                poisECMData_obj$X_list)[, 1]
    }

    if (verbose || (0 == (em_idx %% 100))) {
      cat("beta", beta_tracker[, 1 + em_idx], "\n")
      cat("tau2, gamma",
          cbind(tau2_tracker[, 1 + em_idx], gamma_tracker[, 1 + em_idx]),
          "\n")
      cat("log-lik", data_log_like_tracker[, em_idx], "\n")

      cat(
        paste(
          "Wrapping up of EM Iteration",
          em_idx,
          "of",
          em_iters,
          ";",
          Sys.time() - t0_EM,
          "Elapsed"
        ),
        "\n"
      )
    }

    # Early stopping:
    #    NULL: Don't stop early.
    #    "abs_loglike": Difference in absolute data log likelihood, |old - new|.
    #    "rel_loglike": Relative difference in data log likelihood, |old - new|/|old|.
    #      Data log likelihood taken across alll areas.
    #    "abs_beta_norm": Absolute L2 norm of difference in beta, ||old - new||_2.
    #    "rel_beta_norm": Relative L2 norm of difference in beta,
    #      ||old - new||_2 / ||old||_2.
    if ((em_min_iters < em_idx) & !is.null(em_stopping)) {
      # Likelhood
      ll_old <- sum(data_log_like_tracker[, em_idx - 1], na.rm = TRUE)
      ll_new <- sum(data_log_like_tracker[, em_idx], na.rm = TRUE)
      ll_diff <- (ll_new - ll_old)
      ll_rel_diff <- ll_diff / ll_old

      # Beta
      beta_diff_norm <- sqrt(sum((beta_tracker[, em_idx + 1] - beta_tracker[, em_idx])^2))
      beta_old_norm <- sqrt(sum(beta_tracker[, em_idx]^2))

      # Expected
      e_ll_old <- sum(expected_log_like_tracker[, em_idx - 1], na.rm = TRUE)
      e_ll_new <- sum(expected_log_like_tracker[, em_idx], na.rm = TRUE)
      e_ll_diff <- (e_ll_new - e_ll_old)
      e_ll_rel_diff <- e_ll_diff / e_ll_old
      if (verbose || (0 == (em_idx %% 100))) {
        cat(
          paste0(
            "End of iteration ",
            em_idx,
            "; Norm of difference in beta: ",
            beta_diff_norm,
            ", Relative difference in beta: ",
            beta_diff_norm / beta_old_norm,
            "; Difference in log-likelihood: ",
            ll_diff,
            "; Relative difference in log-likelihood: ",
            ll_rel_diff,
            "; Difference in expected log-likelihood: ",
            e_ll_diff,
            "; Relative difference in expected log-likelihood: ",
            e_ll_rel_diff
          ),
          "\n"
        )
      }

      if ("abs_loglike" == em_stopping) {
        if (abs(ll_diff) < em_tol) {
          cat("Ending early", "\n")
          break
        }
      }
      else if ("rel_loglike" == em_stopping) {
        if (abs(ll_rel_diff) < em_tol) {
          cat("Ending early", "\n")
          break
        }
      }
      else if ("abs_beta_norm" == em_stopping) {
        if (beta_diff_norm < em_tol) {
          cat("Ending early", "\n")
          break
        }
      }
      else if ("rel_beta_norm" == em_stopping) {
        if (beta_diff_norm / beta_old_norm < em_tol) {
          cat("Ending early", "\n")
          break
        }
      }
      else {
        warning("Invalid value for em_stopping.")
      }
    }
  }
  t1_EM = Sys.time()
  cat("Time", t1_EM - t0_EM, "\n")

  # Compute Negative Hessians
  beta_neghessian <- neg_hessian_beta(Q_hat_list, tau2_tracker[, em_idx + 1], poisECMData_obj$X_list)
  tau2_neghessian <- rep(NA, n_areas)
  gamma_neghessian <- rep(NA, n_areas)
  for (area_idx in 1:n_areas) {
    tau2_neghessian[area_idx] <- neg_hessian_tau2(
      Vhat_list[[area_idx]],
      eta_hat_list[[area_idx]],
      Q_hat_list[[area_idx]],
      tau2_tracker[area_idx, em_idx + 1],
      beta_tracker[, em_idx + 1],
      poisECMData_obj$X_list[[area_idx]]
    )

    if ("CAR" == model_type) {
      gamma_neghessian[area_idx] <- neg_hessian_gamma_CAR(gamma_tracker[area_idx, em_idx + 1], eig_val_list[[area_idx]])
    }
    else if ("SAR" == model_type) {
      gamma_neghessian[area_idx] <- neg_hessian_gamma_SAR(
        Vhat_list[[area_idx]],
        eta_hat_list[[area_idx]],
        gamma_tracker[area_idx, em_idx + 1],
        tau2_tracker[area_idx, em_idx + 1],
        beta_tracker[, em_idx + 1],
        poisECMData_obj$X_list[[area_idx]],
        poisECMData_obj$W_list[[area_idx]],
        poisECMData_obj$D_list[[area_idx]],
        eig_val_list[[area_idx]]
      )
    }
    else if ("Leroux" == model_type) {
      gamma_neghessian[area_idx] <- neg_hessian_gamma_Leroux(gamma_tracker[area_idx, em_idx + 1], eig_val_list[[area_idx]])
    }
  }

  ## Name stuff

  rownames(beta_tracker) <- colnames(poisECMData_obj$X_list[[1]])
  rownames(gamma_tracker) <- names(poisECMData_obj$W_list)
  rownames(tau2_tracker) <- names(poisECMData_obj$W_list)

  rownames(R2_tracker) <- names(poisECMData_obj$W_list)
  rownames(MSE_tracker) <- names(poisECMData_obj$W_list)
  rownames(data_log_like_tracker) <- names(poisECMData_obj$W_list)
  rownames(expected_log_like_tracker) <- names(poisECMData_obj$W_list)
  rownames(resid_moran_nb) <- names(poisECMData_obj$W_list)
  dimnames(resid_moran) <- list(
    names(poisECMData_obj$W_list),
    1:dim(resid_moran)[2],
    c("Moran_I", "ExpectedMoran_I", "PValue")
  )

  names(start_idx_list) <- names(poisECMData_obj$W_list)
  names(eig_val_list) <- names(poisECMData_obj$W_list)

  names(tau2_neghessian) <- names(poisECMData_obj$W_list)
  names(gamma_neghessian) <- names(poisECMData_obj$W_list)
  rownames(beta_neghessian) <- colnames(poisECMData_obj$X_list[[1]])
  colnames(beta_neghessian) <- colnames(poisECMData_obj$X_list[[1]])

  rownames(fit_tracker) <- Reduce(c, lapply(poisECMData_obj$counts_list, colnames))
  rownames(eta_tracker) <- rownames(fit_tracker)
  rownames(theta_tracker) <- rownames(fit_tracker)

  # If eigenvalues supplied, don't bother returning them to save space
  eigs_supplied <- (!is.null(poisECMData_obj$eig_CS_list) ||
                      !is.null(poisECMData_obj$eig_L_list))
  if (eigs_supplied) {
    eig_val_list <- list()
  }

  out <- (structure(
    list(
      # Coefficients, spatial parameters
      beta_hat = beta_tracker[, (em_idx + 1)],
      gamma_hat = gamma_tracker[, (em_idx + 1)],
      tau2_hat = tau2_tracker[, (em_idx + 1)],

      # Fitted values and residuals, fitted Poisson parameters, estimated random effects
      predictions = fit_tracker[, em_idx],
      residuals = Reduce(c, z_list) - fit_tracker[, em_idx],
      eta_hat = eta_tracker[, em_idx],
      theta_hat = theta_tracker[, em_idx],
      phi_hat = eta_tracker[, em_idx] - Reduce(rbind, poisECMData_obj$X_list) %*% beta_tracker[, (em_idx + 1)],

      # Full paths of coefficients, spatial parameters
      gamma_tracker = gamma_tracker[, 1:(em_idx + 1)],
      tau2_tracker = tau2_tracker[, 1:(em_idx + 1)],
      beta_tracker = beta_tracker[, 1:(em_idx + 1)],

      # Trackers of various fit diagnostics
      R2_tracker = R2_tracker[, 1:em_idx],
      MSE_tracker = MSE_tracker[, 1:em_idx],
      data_log_like_tracker = data_log_like_tracker[, 1:em_idx],
      expected_log_like_tracker = expected_log_like_tracker[, 1:em_idx],
      resid_moran = resid_moran[, 1:em_idx, ],
      resid_moran_nb = resid_moran_nb[, 1:em_idx],

      # Utilities, just in case we want to reconstruct stuff
      start_idx_list = start_idx_list,
      eig_val_list = eig_val_list,

      # Fit time
      time = difftime(t1_EM, t0_EM),

      # Hessians (for standard errors)
      beta_neghessian = beta_neghessian,
      tau2_neghessian = tau2_neghessian,
      gamma_neghessian = gamma_neghessian,

      run_settings = list(
        gene_name = as.character(gene_name),
        gene_idx = gene_idx,
        model_type = as.character(model_type),
        em_iters = em_iters,
        opt_iters = opt_iters,
        em_min_iters = em_min_iters,
        em_tol = em_tol,
        em_stopping = em_stopping,
        beta_init = beta_init,
        gamma_init = gamma_init,
        tau2_init = tau2_init,
        verbose = verbose,
        em_iters_actual = em_idx,
        eigs_supplied = eigs_supplied
      )
    ),
    class = "poisECMOutput"
  ))

  out$performanceSummary <- summarizePoisECMPerformance(poisECMData_obj, out)
  return(out)
}

#' Check inputs for the poisECM method.
#'  Make sure that the input object has everything needed to run without error.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param poisECMData_obj Object containing data.
#'  Created by the prepData method.
#'
#' @note Does not return anything.
#' @note This method can be used to check a hand-created input object.
#'  E.g., if a user does not want to use prepData.
#'
#' @returns Nothing.
#'
#' @import Matrix
#' @export
checkInputsPoisECM <- function (poisECMData_obj) {
  # Check that the bare minimum is present
  stopifnot(!is.null(poisECMData_obj$counts_list))
  stopifnot(!is.null(poisECMData_obj$W_list))
  stopifnot(!is.null(poisECMData_obj$D_list))
  stopifnot(!is.null(poisECMData_obj$X_list))
  stopifnot(!is.null(poisECMData_obj$library_size_list))

  # Check for optional components
  if (is.null(poisECMData_obj$coords_list)) {
    warning("Coordinates are not present.")
  }
  if (is.null(poisECMData_obj$eig_CS_list) &&
      is.null(poisECMData_obj$eig_L_list)) {
    warning(
      "Eigenvalues are not present; initial computation will be performed when running poisECM (slow)."
    )
  }

  # Counts-specific checks
  # Check that there is at least one gene, really that it's a matrix.
  stopifnot(1 <= min(sapply(poisECMData_obj$counts_list, nrow)))
  # Check that same number of genes present
  stopifnot(1 == length(unique(
    lapply(poisECMData_obj$counts_list, nrow)
  )))
  # Check ordering by rownames (gene names)
  stopifnot(all(sapply(
    lapply(poisECMData_obj$counts_list, rownames),
    identical,
    rownames(poisECMData_obj$counts_list[[1]])
  )))

  # Coordinates-specific checks
  if (!is.null(poisECMData_obj$coords_list)) {
    # Need at least 2 coordinates
    stopifnot(1 < min(sapply(poisECMData_obj$coords_list, ncol)))
  }

  # Check that the number of measurements/cells is the same across lists
  # Implicit check that the number of entries in list is the same (number of samples)
  stopifnot(identical(
    sapply(poisECMData_obj$counts_list, ncol),
    sapply(poisECMData_obj$W_list, nrow)
  ))
  stopifnot(identical(
    sapply(poisECMData_obj$counts_list, ncol),
    sapply(poisECMData_obj$D_list, nrow)
  ))
  stopifnot(identical(
    sapply(poisECMData_obj$counts_list, ncol),
    sapply(poisECMData_obj$X_list, nrow)
  ))
  stopifnot(identical(
    sapply(poisECMData_obj$counts_list, ncol),
    sapply(poisECMData_obj$library_size_list, length)
  ))
  if (!is.null(poisECMData_obj$coords_list)) {
    stopifnot(identical(
      sapply(poisECMData_obj$counts_list, ncol),
      sapply(poisECMData_obj$coords_list, nrow)
    ))
  }
  if (!is.null(poisECMData_obj$eig_CS_list)) {
    stopifnot(identical(
      sapply(poisECMData_obj$counts_list, ncol),
      sapply(poisECMData_obj$eig_CS_list, length)
    ))
  }
  if (!is.null(poisECMData_obj$eig_L_list)) {
    stopifnot(identical(
      sapply(poisECMData_obj$counts_list, ncol),
      sapply(poisECMData_obj$eig_L_list, length)
    ))
  }

  # Check that the names of measurements match
  for (idx in 1:length(poisECMData_obj$counts_list)) {
    stopifnot(identical(
      colnames(poisECMData_obj$counts_list[[idx]]),
      rownames(poisECMData_obj$X_list[[idx]])
    ))
    stopifnot(identical(
      colnames(poisECMData_obj$counts_list[[idx]]),
      names(poisECMData_obj$library_size_list[[idx]])
    ))
    if (!is.null(poisECMData_obj$coords_list)) {
      stopifnot(identical(
        colnames(poisECMData_obj$counts_list[[idx]]),
        rownames(poisECMData_obj$coords_list[[idx]])
      ))
    }
  }

  # Check ordering of lists (sample IDs)
  stopifnot(identical(
    names(poisECMData_obj$counts_list),
    names(poisECMData_obj$W_list)
  ))
  stopifnot(identical(
    names(poisECMData_obj$counts_list),
    names(poisECMData_obj$D_list)
  ))
  stopifnot(identical(
    names(poisECMData_obj$counts_list),
    names(poisECMData_obj$X_list)
  ))
  stopifnot(identical(
    names(poisECMData_obj$counts_list),
    names(poisECMData_obj$library_size_list)
  ))
  if (!is.null(poisECMData_obj$coords_list)) {
    stopifnot(identical(
      names(poisECMData_obj$counts_list),
      names(poisECMData_obj$coords_list)
    ))
  }
  if (!is.null(poisECMData_obj$eig_CS_list)) {
    stopifnot(identical(
      names(poisECMData_obj$counts_list),
      names(poisECMData_obj$eig_CS_list)
    ))
  }
  if (!is.null(poisECMData_obj$eig_L_list)) {
    stopifnot(identical(
      names(poisECMData_obj$counts_list),
      names(poisECMData_obj$eig_L_list)
    ))
  }
}
