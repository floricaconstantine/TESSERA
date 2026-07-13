## Wrapper functions to fit model.
# Dependencies in file: Matrix.
# Dependencies: Functions from utils.R, models.R, E_step.R, M_step.R.
# Dependences from functions in other files not listed: pracma, gstat.
# Rcpp dependencies: calc_moran.cpp.


#' Fit Multi-Sample Poisson Spatial GLMM via ECM Algorithm
#'
#' Fits a multi-sample Poisson spatial generalized linear mixed model (GLMM)
#' using a shared set of fixed effects across all samples while permitting
#' sample-specific spatial random effects (CAR, SAR, or Leroux). Parameter
#' estimation is performed via an Expectation-Conditional Maximization (ECM)
#' algorithm.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param TESSERAData_obj An object containing prepared data, typically
#'   created by \code{\link{prep_data}}.
#' @param gene_name Character: The name of the gene/measurement (row) to fit.
#' @param model_type Character: The spatial model for random effects.
#'   Options are "CAR", "SAR", or "Leroux".
#' @param em_iters Integer: Maximum number of ECM iterations.
#' @param opt_iters Integer: Number of inner CM steps (Conditional Maximization)
#'   per iteration. The algorithm maximizes the expected likelihood for
#'   \eqn{\tau^2}, then \eqn{\gamma} for each area, and finally \eqn{\beta} while
#'   holding other parameters constant.
#' @param em_min_iters Integer: Minimum number of ECM iterations to perform
#'   before allowing early stopping.
#' @param em_tol Numeric: Convergence tolerance for early stopping.
#' @param em_stopping Character: Metric used for early stopping:
#' \itemize{
#'   \item \code{NULL}: No early stopping.
#'   \item "abs_loglike": Absolute change in total log-likelihood.
#'   \item "rel_loglike": Relative change in total log-likelihood.
#'   \item "abs_beta_norm": Absolute \eqn{L_2} norm of the change in \eqn{\beta}.
#'   \item "rel_beta_norm": Relative \eqn{L_2} norm of the change in \eqn{\beta}.
#' }
#' @param beta_init Initial value for \eqn{\beta}. Options: "glm" (fit a Poisson GLM),
#'   "random" (standard normal), or a numeric vector.
#' @param gamma_init Initial value for \eqn{\gamma}. Options: "moran" (absolute
#'   Moran's I of residuals), "random" (standard uniform), or numeric.
#' @param tau2_init Initial value for \eqn{\tau^2}. Options: "lognormal" (approximate
#'   Poisson-lognormal variance), "var" (variance of residuals), "random",
#'   or numeric.
#' @param verbose Logical: Whether to print iteration-wise parameter updates.
#' @param dense_matrices Logical: If \code{TRUE}, treats the precision matrix \eqn{Q}
#'   as dense during specific E-step calculations. This increases memory usage
#'   but can improve computation speed by 10 to 20 percent.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \strong{beta_hat}: Estimated fixed effect coefficients.
#'   \item \strong{gamma_hat}: Estimated spatial correlation parameters (per area).
#'   \item \strong{tau2_hat}: Estimated spatial scale parameters (per area).
#'   \item \strong{phi_hat}: Estimated spatial random effects.
#'   \item \strong{theta_hat}: Estimated Poisson rate parameters \eqn{exp(X\beta + \phi)}.
#'   \item \strong{beta_tracker}: Matrix of \eqn{\beta} estimates across iterations.
#'   \item \strong{gamma_tracker}: Matrix of \eqn{\gamma} estimates across iterations.
#'   \item \strong{tau2_tracker}: Matrix of \eqn{\tau^2} estimates across iterations.
#'   \item \strong{performanceSummary}: A summary data frame for each sample.
#'   \item \strong{time}: Total execution time.
#' }
#'
#' @references Meng, Xiao-Li, and Donald B. Rubin. "Maximum likelihood estimation via the ECM algorithm: A general framework." Biometrika 80.2 (1993): 267-278.
#'
#' @note For CAR and SAR models, ensure the adjacency structure \eqn{W} contains no
#'   isolated points (every observation must have \eqn{\ge 1} neighbor) to ensure
#'   invertibility.
#'
#' @import Matrix
#' @importFrom Matrix diag
#' @importFrom stats coef cor dpois poisson predict rnorm runif var
#' @importFrom Rcpp sourceCpp evalCpp
#' @importFrom methods as
#' @useDynLib TESSERA
#' @export
#'
#' @examples
#' # Locate the prepped TESSERAData object in inst/extdata
#' rds_path <- system.file("extdata", "example_prepData.rds", package = "TESSERA")
#' # Load the TESSERAData object
#' TESSERA_data <- readRDS(rds_path)
#'
#' # Fit the Poisson generalized spatial linear model
#' TESSERA_out_Leroux <- suppressMessages(suppressWarnings(
#'   TESSERA_lattice(
#'     TESSERAData_obj = TESSERA_data,
#'     gene_name = "example",
#'     model_type = "Leroux",
#'     em_iters = 2,
#'     opt_iters = 1,
#'     verbose = FALSE
#'   )
#' ))
TESSERA_lattice <- function(TESSERAData_obj,
                            gene_name,
                            model_type = "Leroux",
                            em_iters = 200,
                            opt_iters = 5,
                            em_min_iters = 15,
                            em_tol = 1e-3,
                            em_stopping = NULL,
                            beta_init = "glm",
                            gamma_init = "moran",
                            tau2_init = "var",
                            verbose = TRUE,
                            dense_matrices = FALSE) {
  # Start clock
  t0_EM <- Sys.time()
  
  # Check inputs
  check_inputs_TESSERA(TESSERAData_obj)
  
  # Extract gene of interest and associated counts
  gene_idx <- which(rownames(TESSERAData_obj$counts_list[[1]]) == gene_name)
  z_list <- lapply(TESSERAData_obj$counts_list, function (x) {
    x[gene_idx, ]
  })
  
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
  
  # gamma optimization function
  if ("CAR" == model_type) {
    gamma_M_fcn <- M_step_gamma_CAR
  } else if ("SAR" == model_type) {
    gamma_M_fcn <- M_step_gamma_SAR
  } else if ("Leroux" == model_type) {
    gamma_M_fcn <- M_step_gamma_Leroux
  } else {
    stop("Invalid model_type")
  }
  
  # Number of areas
  n_areas <- length(TESSERAData_obj$W_list)
  # Total number of points
  n_total_points <- 0
  # Starting points in larger matrices
  start_idx_list <- rep(NA, length(z_list))
  for (idx in 1:length(TESSERAData_obj$X_list)) {
    n_total_points <- n_total_points + nrow(TESSERAData_obj$X_list[[idx]])
    if (1 == idx) {
      start_idx_list[idx] <- 1
    } else {
      start_idx_list[idx] <- start_idx_list[idx - 1] + length(z_list[[idx - 1]])
    }
  }
  # Dimension of coefficient vector
  beta_dim <- ncol(TESSERAData_obj$X_list[[1]])
  
  # Identify which specific eigenvalues we actually need
  needs_CS <- (model_type %in% c("CAR", "SAR"))
  needs_L  <- (model_type == "Leroux")
  
  # Trigger computation if the SPECIFIC list required is missing
  if ((needs_CS && is.null(TESSERAData_obj$eig_CS_list)) ||
      (needs_L && is.null(TESSERAData_obj$eig_L_list))) {
    if (verbose)
      message("Required eigenvalues not supplied: Computing now.")
    eig_val_list <- list()
    
    for (area_idx in 1:n_areas) {
      if (needs_CS) {
        eig_val_list[[area_idx]] <- Re(eigen(
          Matrix::solve(TESSERAData_obj$D_list[[area_idx]], TESSERAData_obj$W_list[[area_idx]]),
          FALSE,
          only.values = TRUE
        )$values)
      } else {
        eig_val_list[[area_idx]] <- Re(
          eigen(
            TESSERAData_obj$D_list[[area_idx]] - TESSERAData_obj$W_list[[area_idx]],
            TRUE,
            only.values = TRUE
          )$values
        )
      }
    }
  } else {
    # If they already exist, grab the correct set
    eig_val_list <- if (needs_L)
      TESSERAData_obj$eig_L_list
    else
      TESSERAData_obj$eig_CS_list
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
    message("Random initialization for beta.", "\n")
  } else if (is.character(beta_init) && ("glm" == beta_init)) {
    # Let's be intelligent: initialize with a basic GLM
    
    # Stack into a single vector/matrix
    # unlist(..., use.names=FALSE) is extremely fast and avoids metadata overhead
    z_vec   <- unlist(z_list, use.names = FALSE)
    lib_vec <- unlist(TESSERAData_obj$library_size_list, use.names = FALSE)
    # do.call(rbind, ...) is the standard way to merge a list of matrices
    # with minimal memory fragmentation.
    X_mat   <- do.call(rbind, TESSERAData_obj$X_list)
    
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
    
    message("GLM initialization for beta.", "\n")
  } else if (is.numeric(beta_init) &&
             (is.vector(beta_init) || is.matrix(beta_init))) {
    # Pass in a value
    beta_tracker[, 1] <- as.vector(beta_init)
    
    message("Pre-defined initialization for beta.", "\n")
  } else {
    stop("Invalid initialization for beta.")
  }
  # Handle NA
  beta_tracker[is.nan(beta_tracker[, 1]), 1] <- 0.0
  beta_tracker[is.na(beta_tracker[, 1]), 1] <- 0.0
  message("Initial beta ", paste(beta_tracker[, 1], collapse = " "), "\n")
  
  # Initialize parameters: gamma
  if (is.character(gamma_init) && ("moran" == gamma_init)) {
    for (area_idx in 1:n_areas) {
      tmp_z <- log((0.5 + z_list[[area_idx]]) / TESSERAData_obj$library_size_list[[area_idx]]) - TESSERAData_obj$X_list[[area_idx]] %*% beta_tracker[, 1]
      gamma_tracker[area_idx, 1] <- abs(moran_I_nb(tmp_z, TESSERAData_obj$W_list[[area_idx]]))
    }
    message("Moran initialization for gamma", "\n")
  } else if (is.character(gamma_init) && ("random" == gamma_init)) {
    gamma_tracker[, 1] <- stats::runif(n_areas)
    message("Random initialization for gamma", "\n")
  } else if (is.numeric(gamma_init)) {
    gamma_tracker[, 1] <- gamma_init
    message("Predefined initialization for gamma", "\n")
  } else {
    stop("Invalid initialization for gamma.")
  }
  # Handle zero counts/NaN values
  gamma_tracker[is.nan(gamma_tracker[, 1]), 1] <- 0.0
  gamma_tracker[is.na(gamma_tracker[, 1]), 1] <- 0.0
  
  message("Initial gamma ", paste(gamma_tracker[, 1], collapse = " "), "\n")
  
  
  
  # Initialize parameters: tau^2
  if (is.character(tau2_init) && ("lognormal" == tau2_init)) {
    for (area_idx in 1:n_areas) {
      glm_out <- stats::glm((z_list[[area_idx]] / TESSERAData_obj$library_size_list[[area_idx]]) ~ 0 + TESSERAData_obj$X_list[[area_idx]],
                            family = stats::poisson()
      )
      tau2_tracker[area_idx, 1] <- max(
        1e-6,
        2 * log(mean(
          (z_list[[area_idx]] / TESSERAData_obj$library_size_list[[area_idx]])
          / stats::predict(glm_out, type = "response")
        )) #- 2
      )
    }
    message("Log-Normal initialization for tau^2", "\n")
  } else if (is.character(tau2_init) && ("var" == tau2_init)) {
    for (area_idx in 1:n_areas) {
      tmp_z <- log((0.5 + z_list[[area_idx]]) / TESSERAData_obj$library_size_list[[area_idx]]) - TESSERAData_obj$X_list[[area_idx]] %*% beta_tracker[, 1]
      tau2_tracker[area_idx, 1] <- stats::var(tmp_z)
    }
    message("Variance initialization for tau^2", "\n")
  } else if (is.character(tau2_init) && ("random" == tau2_init)) {
    tau2_tracker[, 1] <- abs(stats::rnorm(n_areas))
    message("Random initialization for tau^2", "\n")
  } else if (is.numeric(tau2_init)) {
    tau2_tracker[, 1] <- tau2_init
    message("Predefined initialization for tau^2", "\n")
  } else {
    stop("Invalid initialization for tau^2.")
  }
  # Handle zero counts/NaN values
  tau2_tracker[is.nan(tau2_tracker[, 1]), 1] <- 0.0
  tau2_tracker[is.na(tau2_tracker[, 1]), 1] <- 0.0
  tau2_tracker[tau2_tracker[, 1] == 0.0, 1] <- 1e-3
  message("Initial tau2 ", paste(tau2_tracker[, 1], collapse = " "), "\n")
  
  # Also initialize dependency Q as a function of W, D, and gamma
  Q_hat_list <- list()
  for (area_idx in 1:n_areas) {
    Q_hat_list[[area_idx]] <- Q_fcn(TESSERAData_obj$W_list[[area_idx]],
                                    TESSERAData_obj$D_list[[area_idx]],
                                    gamma_tracker[area_idx, 1])
  }
  
  # Loop over EM iterations
  for (em_idx in 1:em_iters) {
    if (verbose || (0 == (em_idx %% 100))) {
      message(
        paste(
          "Start of EM Iteration ",
          em_idx,
          " of ",
          em_iters,
          "; ",
          Sys.time() - t0_EM,
          " Elapsed"
        ),
        "\n"
      )
    }
    
    
    # Run E-Step: Get covariance and mean of eta
    trace_scalars_list <- list() # REPLACES Vhat_list
    eta_hat_list <- list()
    
    for (area_idx in 1:n_areas) {
      # Compute Vhat for the CURRENT area only
      if (dense_matrices) {
        Vhat_current <- E_step_Vhat(as.matrix(Q_hat_list[[area_idx]]), tau2_tracker[area_idx, em_idx], z_list[[area_idx]])
      } else {
        Vhat_current <- E_step_Vhat(Q_hat_list[[area_idx]], tau2_tracker[area_idx, em_idx], z_list[[area_idx]])
      }
      
      # Get the mean (Vhat argument removed; reconstructs sparse Vinv internally)
      eta_hat_list[[area_idx]] <- E_step_etahat(
        Q_hat_list[[area_idx]],
        tau2_tracker[area_idx, em_idx],
        beta_tracker[, em_idx],
        TESSERAData_obj$X_list[[area_idx]],
        z_list[[area_idx]],
        TESSERAData_obj$library_size_list[[area_idx]]
      )
      
      # Get a few more things out of the E-step: theta and predictions
      theta_hat <- as.numeric(E_step_thetahat(Vhat_current, eta_hat_list[[area_idx]]))
      z_hat <- as.numeric(E_step_predict(theta_hat, TESSERAData_obj$library_size_list[[area_idx]]))
      
      # Store stuff
      fit_tracker[start_idx_list[area_idx]:(start_idx_list[area_idx] + length(z_hat) - 1), em_idx] <- z_hat
      eta_tracker[start_idx_list[area_idx]:(start_idx_list[area_idx] + length(eta_hat_list[[area_idx]]) - 1), em_idx] <- as.numeric(eta_hat_list[[area_idx]])
      theta_tracker[start_idx_list[area_idx]:(start_idx_list[area_idx] + length(theta_hat) - 1), em_idx] <- theta_hat
      
      # Performance
      R2_tracker[area_idx, em_idx] <- stats::cor(z_hat, z_list[[area_idx]])^2
      MSE_tracker[area_idx, em_idx] <- mean(abs(z_hat - z_list[[area_idx]])^2)
      
      # Observed, expected, sd
      if (!(
        is.null(TESSERAData_obj$coords_list[[area_idx]][, 1]) ||
        is.null(TESSERAData_obj$coords_list[[area_idx]][, 2])
      )) {
        resid_moran[area_idx, em_idx, ] <- calc_moran(
          z_hat - z_list[[area_idx]],
          TESSERAData_obj$coords_list[[area_idx]][, 1],
          TESSERAData_obj$coords_list[[area_idx]][, 2]
        )
      }
      resid_moran_nb[area_idx, em_idx] <- moran_I_nb(z_hat - z_list[[area_idx]], TESSERAData_obj$W_list[[area_idx]])
      
      # Data log likelihood
      data_log_like_tracker[area_idx, em_idx] <- sum(
        stats::dpois(
          round(z_list[[area_idx]]),
          theta_hat * TESSERAData_obj$library_size_list[[area_idx]],
          log = TRUE
        ),
        na.rm = TRUE
      )
      
      # Expected log likelihood (Takes Vhat_current before it gets deleted!)
      expected_log_like_tracker[area_idx, em_idx] <- expected_loglike(
        Vhat_current,
        eta_hat_list[[area_idx]],
        Q_hat_list[[area_idx]],
        gamma_tracker[area_idx, em_idx],
        tau2_tracker[area_idx, em_idx],
        beta_tracker[, em_idx],
        TESSERAData_obj$X_list[[area_idx]],
        TESSERAData_obj$W_list[[area_idx]],
        TESSERAData_obj$D_list[[area_idx]],
        eig_val_list[[area_idx]],
        model_type
      )
      
      # Extract Trace Scalars for the M-Step
      W_mat <- TESSERAData_obj$W_list[[area_idx]]
      D_mat <- TESSERAData_obj$D_list[[area_idx]]
      
      if ("CAR" == model_type) {
        # Old: W_sparse <- methods::as(W_mat, "dgCMatrix")
        W_sparse <- methods::as(methods::as(W_mat, "generalMatrix"), "CsparseMatrix")
        idx_W <- cbind(W_sparse@i + 1L, rep(1:ncol(W_sparse), diff(W_sparse@p)))
        
        trace_scalars_list[[area_idx]] <- list(
          tr_DV = sum(Matrix::diag(D_mat) * Matrix::diag(Vhat_current)),
          tr_WV = sum(W_sparse@x * Vhat_current[idx_W])
        )
        
      } else if ("SAR" == model_type) {
        # Old: W_sparse <- methods::as(W_mat, "dgCMatrix")
        W_sparse <- methods::as(methods::as(W_mat, "generalMatrix"), "CsparseMatrix")
        idx_W <- cbind(W_sparse@i + 1L, rep(1:ncol(W_sparse), diff(W_sparse@p)))
        
        D_inv <- Matrix::Diagonal(dim(W_mat)[1], 1 / Matrix::diag(D_mat))
        # Old: WZ_sparse <- methods::as(W_mat %*% (D_inv %*% W_mat), "dgCMatrix")
        WZ_sparse <- methods::as(methods::as(W_mat %*% (D_inv %*% W_mat), "generalMatrix"), "CsparseMatrix")
        idx_WZ <- cbind(WZ_sparse@i + 1L, rep(1:ncol(WZ_sparse), diff(WZ_sparse@p)))
        
        trace_scalars_list[[area_idx]] <- list(
          tr_DV = sum(Matrix::diag(D_mat) * Matrix::diag(Vhat_current)),
          tr_WV = sum(W_sparse@x * Vhat_current[idx_W]),
          tr_WZV = sum(WZ_sparse@x * Vhat_current[idx_WZ])
        )
        
      } else if ("Leroux" == model_type) {
        # Old: DWI_sparse <- methods::as(D_mat - W_mat - Matrix::Diagonal(...), "dgCMatrix")
        DWI_sparse <- methods::as(methods::as(D_mat - W_mat - Matrix::Diagonal(dim(W_mat)[1], 1), "generalMatrix"), "CsparseMatrix")
        idx_DWI <- cbind(DWI_sparse@i + 1L, rep(1:ncol(DWI_sparse), diff(DWI_sparse@p)))
        
        trace_scalars_list[[area_idx]] <- list(
          tr_V = sum(Matrix::diag(Vhat_current)),
          tr_DWIV = sum(DWI_sparse@x * Vhat_current[idx_DWI])
        )
      }
      
      # DESTROY the sparse Vhat to prevent Memory Overflow
      rm(Vhat_current)
      gc() # Force R to collect the garbage immediately
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
          trace_scalars_list[[area_idx]],
          gamma_tracker[area_idx, 1 + em_idx],
          # Current gamma
          model_type,
          eta_hat_list[[area_idx]],
          Q_hat_list[[area_idx]],
          beta_tracker[, 1 + em_idx],
          TESSERAData_obj$X_list[[area_idx]]
        )
        
        # Optimize in gamma
        gamma_out <- gamma_M_fcn(
          trace_scalars_list[[area_idx]],
          # Replaced Vhat_list
          eta_hat_list[[area_idx]],
          tau2_tracker[area_idx, 1 + em_idx],
          beta_tracker[, 1 + em_idx],
          TESSERAData_obj$X_list[[area_idx]],
          TESSERAData_obj$W_list[[area_idx]],
          TESSERAData_obj$D_list[[area_idx]],
          eig_val_list[[area_idx]],
          gamma_tracker[area_idx, 1 + em_idx]
        )
        
        # Handle optimization going off the rails
        gamma_val <- gamma_out$gamma_hat # Use local variable to avoid masking
        
        if (("CAR" == model_type) | ("SAR" == model_type)) {
          if (-1 > gamma_val) {
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
            gamma_val <- max(-1, gamma_val)
          }
          if (1 < gamma_val) {
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
            gamma_val <- min(1, gamma_val)
          }
        } else if ("Leroux" == model_type) {
          if (0 > gamma_val) {
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
            gamma_val <- max(0, gamma_val)
          }
          if (1 < gamma_val) {
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
            gamma_val <- min(1, gamma_val)
          }
        }
        gamma_tracker[area_idx, 1 + em_idx] <- gamma_val
        
        # Update precision matrix Q
        Q_hat_list[[area_idx]] <- Q_fcn(TESSERAData_obj$W_list[[area_idx]],
                                        TESSERAData_obj$D_list[[area_idx]],
                                        gamma_tracker[area_idx, em_idx + 1])
      }
      
      # Optimize in beta
      beta_tracker[, 1 + em_idx] <- M_step_beta(eta_hat_list,
                                                Q_hat_list,
                                                tau2_tracker[, 1 + em_idx],
                                                TESSERAData_obj$X_list)[, 1]
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
        message(
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
          message("Ending early", "\n")
          break
        }
      }
      else if ("rel_loglike" == em_stopping) {
        if (abs(ll_rel_diff) < em_tol) {
          message("Ending early", "\n")
          break
        }
      }
      else if ("abs_beta_norm" == em_stopping) {
        if (beta_diff_norm < em_tol) {
          message("Ending early", "\n")
          break
        }
      }
      else if ("rel_beta_norm" == em_stopping) {
        if (beta_diff_norm / beta_old_norm < em_tol) {
          message("Ending early", "\n")
          break
        }
      }
      else {
        warning("Invalid value for em_stopping.")
      }
    }
  }
  t1_EM <- Sys.time()
  message("Time ", t1_EM - t0_EM, "\n")
  
  # Compute Negative Hessians
  beta_neghessian <- neg_hessian_beta(Q_hat_list, tau2_tracker[, em_idx + 1], TESSERAData_obj$X_list)
  tau2_neghessian <- rep(NA, n_areas)
  gamma_neghessian <- rep(NA, n_areas)
  
  for (area_idx in 1:n_areas) {
    # Recompute Vhat for the final parameter estimates to calculate Hessians
    if (dense_matrices) {
      Vhat_current <- E_step_Vhat(as.matrix(Q_hat_list[[area_idx]]), tau2_tracker[area_idx, em_idx + 1], z_list[[area_idx]])
    } else {
      Vhat_current <- E_step_Vhat(Q_hat_list[[area_idx]], tau2_tracker[area_idx, em_idx + 1], z_list[[area_idx]])
    }
    
    # RECOMPUTE eta_hat for the final parameter estimates
    eta_hat_final <- E_step_etahat(
      Q_hat_list[[area_idx]],
      tau2_tracker[area_idx, em_idx + 1],
      beta_tracker[, em_idx + 1],
      TESSERAData_obj$X_list[[area_idx]],
      z_list[[area_idx]],
      TESSERAData_obj$library_size_list[[area_idx]]
    )
    
    # Compute tau2 negative Hessian
    tau2_neghessian[area_idx] <- neg_hessian_tau2(
      Vhat_current,
      eta_hat_final,
      Q_hat_list[[area_idx]],
      tau2_tracker[area_idx, em_idx + 1],
      beta_tracker[, em_idx + 1],
      TESSERAData_obj$X_list[[area_idx]]
    )
    
    # Compute gamma negative Hessian
    if ("CAR" == model_type) {
      gamma_neghessian[area_idx] <- neg_hessian_gamma_CAR(gamma_tracker[area_idx, em_idx + 1], eig_val_list[[area_idx]])
    } else if ("SAR" == model_type) {
      gamma_neghessian[area_idx] <- neg_hessian_gamma_SAR(
        Vhat_current,
        eta_hat_final,
        gamma_tracker[area_idx, em_idx + 1],
        tau2_tracker[area_idx, em_idx + 1],
        beta_tracker[, em_idx + 1],
        TESSERAData_obj$X_list[[area_idx]],
        TESSERAData_obj$W_list[[area_idx]],
        TESSERAData_obj$D_list[[area_idx]],
        eig_val_list[[area_idx]]
      )
    } else if ("Leroux" == model_type) {
      gamma_neghessian[area_idx] <- neg_hessian_gamma_Leroux(gamma_tracker[area_idx, em_idx + 1], eig_val_list[[area_idx]])
    }
    
    # Destroy Vhat to maintain low memory footprint
    rm(Vhat_current)
    gc()
  }
  
  
  ## Name stuff
  
  rownames(beta_tracker) <- colnames(TESSERAData_obj$X_list[[1]])
  rownames(gamma_tracker) <- names(TESSERAData_obj$W_list)
  rownames(tau2_tracker) <- names(TESSERAData_obj$W_list)
  
  rownames(R2_tracker) <- names(TESSERAData_obj$W_list)
  rownames(MSE_tracker) <- names(TESSERAData_obj$W_list)
  rownames(data_log_like_tracker) <- names(TESSERAData_obj$W_list)
  rownames(expected_log_like_tracker) <- names(TESSERAData_obj$W_list)
  rownames(resid_moran_nb) <- names(TESSERAData_obj$W_list)
  dimnames(resid_moran) <- list(
    names(TESSERAData_obj$W_list),
    1:dim(resid_moran)[2],
    c("Moran_I", "ExpectedMoran_I", "PValue")
  )
  
  names(start_idx_list) <- names(TESSERAData_obj$W_list)
  names(eig_val_list) <- names(TESSERAData_obj$W_list)
  
  names(tau2_neghessian) <- names(TESSERAData_obj$W_list)
  names(gamma_neghessian) <- names(TESSERAData_obj$W_list)
  rownames(beta_neghessian) <- colnames(TESSERAData_obj$X_list[[1]])
  colnames(beta_neghessian) <- colnames(TESSERAData_obj$X_list[[1]])
  
  rownames(fit_tracker) <- Reduce(c, lapply(TESSERAData_obj$counts_list, colnames))
  rownames(eta_tracker) <- rownames(fit_tracker)
  rownames(theta_tracker) <- rownames(fit_tracker)
  
  # If eigenvalues supplied, don't bother returning them to save space
  eigs_supplied <- (!is.null(TESSERAData_obj$eig_CS_list) ||
                      !is.null(TESSERAData_obj$eig_L_list))
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
      phi_hat = eta_tracker[, em_idx] - do.call(rbind, TESSERAData_obj$X_list) %*% beta_tracker[, (em_idx + 1)],
      
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
    class = "TESSERAOutput"
  ))
  
  out$performanceSummary <- summarize_TESSERA(TESSERAData_obj, out)
  return(out)
}

#' Check inputs for the TESSERA method.
#'  Make sure that the input object has everything needed to run without error.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param TESSERAData_obj Object containing data.
#'  Created by the prep_data method.
#'
#' @note Does not return anything.
#' @note This method can be used to check a hand-created input object.
#'  E.g., if a user does not want to use prep_data.
#'
#' @returns Nothing.
#'
#' @import Matrix
check_inputs_TESSERA <- function (TESSERAData_obj) {
  # Check that the bare minimum is present
  stopifnot(!is.null(TESSERAData_obj$counts_list))
  stopifnot(!is.null(TESSERAData_obj$W_list))
  stopifnot(!is.null(TESSERAData_obj$D_list))
  stopifnot(!is.null(TESSERAData_obj$X_list))
  stopifnot(!is.null(TESSERAData_obj$library_size_list))
  
  # Check for optional components
  if (is.null(TESSERAData_obj$coords_list)) {
    warning("Coordinates are not present.")
  }
  if (is.null(TESSERAData_obj$eig_CS_list) &&
      is.null(TESSERAData_obj$eig_L_list)) {
    warning(
      "Eigenvalues are not present; initial computation will be performed when running TESSERA (slow)."
    )
  }
  
  # Counts-specific checks
  # Check that there is at least one gene, really that it's a matrix.
  stopifnot(1 <= min(sapply(TESSERAData_obj$counts_list, nrow)))
  # Check that same number of genes present
  stopifnot(1 == length(unique(
    lapply(TESSERAData_obj$counts_list, nrow)
  )))
  # Check ordering by rownames (gene names)
  stopifnot(all(sapply(
    lapply(TESSERAData_obj$counts_list, rownames),
    identical,
    rownames(TESSERAData_obj$counts_list[[1]])
  )))
  
  # Coordinates-specific checks
  if (!is.null(TESSERAData_obj$coords_list)) {
    # Need at least 2 coordinates
    stopifnot(1 < min(sapply(TESSERAData_obj$coords_list, ncol)))
  }
  
  # Check that the number of measurements/cells is the same across lists
  # Implicit check that the number of entries in list is the same (number of samples)
  stopifnot(identical(
    sapply(TESSERAData_obj$counts_list, ncol),
    sapply(TESSERAData_obj$W_list, nrow)
  ))
  stopifnot(identical(
    sapply(TESSERAData_obj$counts_list, ncol),
    sapply(TESSERAData_obj$D_list, nrow)
  ))
  stopifnot(identical(
    sapply(TESSERAData_obj$counts_list, ncol),
    sapply(TESSERAData_obj$X_list, nrow)
  ))
  stopifnot(identical(
    sapply(TESSERAData_obj$counts_list, ncol),
    sapply(TESSERAData_obj$library_size_list, length)
  ))
  if (!is.null(TESSERAData_obj$coords_list)) {
    stopifnot(identical(
      sapply(TESSERAData_obj$counts_list, ncol),
      sapply(TESSERAData_obj$coords_list, nrow)
    ))
  }
  if (!is.null(TESSERAData_obj$eig_CS_list)) {
    stopifnot(identical(
      sapply(TESSERAData_obj$counts_list, ncol),
      sapply(TESSERAData_obj$eig_CS_list, length)
    ))
  }
  if (!is.null(TESSERAData_obj$eig_L_list)) {
    stopifnot(identical(
      sapply(TESSERAData_obj$counts_list, ncol),
      sapply(TESSERAData_obj$eig_L_list, length)
    ))
  }
  
  # Check that the names of measurements match
  for (idx in 1:length(TESSERAData_obj$counts_list)) {
    stopifnot(identical(
      colnames(TESSERAData_obj$counts_list[[idx]]),
      rownames(TESSERAData_obj$X_list[[idx]])
    ))
    stopifnot(identical(
      colnames(TESSERAData_obj$counts_list[[idx]]),
      names(TESSERAData_obj$library_size_list[[idx]])
    ))
    if (!is.null(TESSERAData_obj$coords_list)) {
      stopifnot(identical(
        colnames(TESSERAData_obj$counts_list[[idx]]),
        rownames(TESSERAData_obj$coords_list[[idx]])
      ))
    }
  }
  
  # Check ordering of lists (sample IDs)
  stopifnot(identical(
    names(TESSERAData_obj$counts_list),
    names(TESSERAData_obj$W_list)
  ))
  stopifnot(identical(
    names(TESSERAData_obj$counts_list),
    names(TESSERAData_obj$D_list)
  ))
  stopifnot(identical(
    names(TESSERAData_obj$counts_list),
    names(TESSERAData_obj$X_list)
  ))
  stopifnot(identical(
    names(TESSERAData_obj$counts_list),
    names(TESSERAData_obj$library_size_list)
  ))
  if (!is.null(TESSERAData_obj$coords_list)) {
    stopifnot(identical(
      names(TESSERAData_obj$counts_list),
      names(TESSERAData_obj$coords_list)
    ))
  }
  if (!is.null(TESSERAData_obj$eig_CS_list)) {
    stopifnot(identical(
      names(TESSERAData_obj$counts_list),
      names(TESSERAData_obj$eig_CS_list)
    ))
  }
  if (!is.null(TESSERAData_obj$eig_L_list)) {
    stopifnot(identical(
      names(TESSERAData_obj$counts_list),
      names(TESSERAData_obj$eig_L_list)
    ))
  }
}
