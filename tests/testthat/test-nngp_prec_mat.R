library(testthat)
library(Matrix)

# --- Shared Setup ---
set.seed(2026)
n_pts <- 30
nngp_k <- 5
coords <- cbind(runif(n_pts), runif(n_pts))
coords <- coords[order(coords[, 1], coords[, 2]), , drop = FALSE]
sp_dist <- sparse_dist(coords, nngp_k)

# Reference R implementation of neighbor distance list
build_nb_dist_ref <- function(coords, sp_dist, nngp_k) {
  n_cols <- ncol(sp_dist)
  area_nb_dists <- vector("list", n_cols)
  for (idx in 1:n_cols) {
    keep_idx <- which(!is.na(sp_dist[(1 + nngp_k):(2 * nngp_k), idx]))
    if (length(keep_idx) > 1) {
      area_nb_dists[[idx]] <- unname(as.matrix(stats::dist(
        coords[sp_dist[nngp_k + keep_idx, idx], , drop = FALSE], diag = TRUE, upper = TRUE
      )))
    } else {
      area_nb_dists[[idx]] <- matrix(0.0, 1, 1)
    }
  }
  area_nb_dists
}

# Reference R implementation of nngp_prec_mat
nngp_prec_mat_ref <- function(sp_dist, nb_dist, cov_type, cov_params) {
  nngp_k <- nrow(sp_dist) / 2
  n_cols <- ncol(sp_dist)
  A_col <- vector("list", n_cols)
  A_row <- vector("list", n_cols)
  A_val <- vector("list", n_cols)
  Dinv_val <- rep(0, n_cols + 1)
  
  for (idx in 1:n_cols) {
    keep_idx <- which(!is.na(sp_dist[(1 + nngp_k):(2 * nngp_k), idx]))
    n_keep <- length(keep_idx)
    if (0 == n_keep)
      next
    
    cov_nb <- sp_dist[keep_idx, idx]
    cov_nb_mat <- nb_dist[[idx]]
    
    if ("Mat" == cov_type) {
      cov_nb <- kernel_matern(cov_nb, cov_params[2], cov_params[3], cov_params[4])
      cov_nb_mat <- kernel_matern(cov_nb_mat, cov_params[2], cov_params[3], cov_params[4])
    } else if ("Exp" == cov_type) {
      cov_nb <- kernel_exp(cov_nb, cov_params[2], cov_params[3])
      cov_nb_mat <- kernel_exp(cov_nb_mat, cov_params[2], cov_params[3])
    } else if ("Gau" == cov_type) {
      cov_nb <- kernel_gauss(cov_nb, cov_params[2], cov_params[3])
      cov_nb_mat <- kernel_gauss(cov_nb_mat, cov_params[2], cov_params[3])
    } else if ("Sph" == cov_type) {
      cov_nb <- kernel_sph(cov_nb, cov_params[2], cov_params[3])
      cov_nb_mat <- kernel_sph(cov_nb_mat, cov_params[2], cov_params[3])
    }
    
    diag(cov_nb_mat) <- diag(cov_nb_mat) + cov_params[1]
    vector_term <- solve(cov_nb_mat, cov_nb)
    
    A_row[[idx]] <- rep(idx + 1, n_keep)
    A_col[[idx]] <- sp_dist[nngp_k + keep_idx, idx]
    A_val[[idx]] <- as.vector(vector_term)
    Dinv_val[idx] <- 1.0 / (cov_params[1] + cov_params[2] - sum(cov_nb * vector_term))
  }
  Dinv_val[length(Dinv_val)] <- 1.0 / (cov_params[1] + cov_params[2])
  
  A_mat <- Matrix::sparseMatrix(
    i = unlist(A_row),
    j = unlist(A_col),
    x = unlist(A_val),
    dims = 1 + c(n_cols, n_cols),
    triangular = TRUE
  )
  id_mat <- Matrix::Diagonal(dim(A_mat)[1], 1)
  I_minus_A <- id_mat - A_mat
  Q <- Matrix::crossprod(I_minus_A, Matrix::Diagonal(length(Dinv_val), Dinv_val) %*% I_minus_A)
  list(Q = Q, Dinv = Dinv_val, A = A_mat)
}

# --- Tests ---

test_that("build_nb_dist_cpp matches base R stats::dist loop exactly", {
  nb_cpp <- build_nb_dist_cpp(coords, sp_dist, nngp_k)
  nb_ref <- build_nb_dist_ref(coords, sp_dist, nngp_k)
  
  expect_length(nb_cpp, ncol(sp_dist))
  expect_equal(nb_cpp, nb_ref, tolerance = 1e-14)
  
  # Bounds checking on malformed inputs
  expect_error(build_nb_dist_cpp(coords, sp_dist, nngp_k = 100),
               "fewer than 2 \\* nngp_k rows")
})

test_that("nngp_prec_triplets_cpp and nngp_prec_mat match reference R across all kernels",
          {
            nb_dist <- build_nb_dist_cpp(coords, sp_dist, nngp_k)
            
            kernels <- list(
              Exp = c(0.1, 1.2, 0.3),
              Gau = c(0.2, 0.8, 0.25),
              Sph = c(0.15, 1.0, 0.4),
              Mat = c(0.1, 1.1, 0.3, 1.5)
            )
            
            for (k_name in names(kernels)) {
              params <- kernels[[k_name]]
              res_cpp <- nngp_prec_mat(sp_dist, nb_dist, k_name, params)
              res_ref <- nngp_prec_mat_ref(sp_dist, nb_dist, k_name, params)
              
              expect_equal(res_cpp$Dinv, res_ref$Dinv, tolerance = 1e-12)
              expect_equal(as.matrix(res_cpp$A), as.matrix(res_ref$A), tolerance = 1e-12)
              expect_equal(as.matrix(res_cpp$Q), as.matrix(res_ref$Q), tolerance = 1e-12)
              expect_true(all(res_cpp$Dinv > 0))
            }
            
            # Error handling on invalid kernel or short parameter vector
            expect_error(
              nngp_prec_triplets_cpp(sp_dist, nb_dist, "InvalidKernel", c(0.1, 1.0, 0.3)),
              "Invalid or unimplemented covariance structure"
            )
            expect_error(
              nngp_prec_triplets_cpp(sp_dist, nb_dist, "Exp", c(0.1, 1.0)),
              "cov_params must have length >= 3"
            )
          })

test_that("fast_trace_quad_cpp matches exact sparse trace Tr(Q * Vhat)", {
  nb_dist <- build_nb_dist_cpp(coords, sp_dist, nngp_k)
  prec <- nngp_prec_mat(sp_dist, nb_dist, "Exp", c(0.1, 1.0, 0.3))
  
  Q <- prec$Q
  Vinv <- Q + Matrix::Diagonal(nrow(Q), runif(nrow(Q), 0.5, 5.0))
  Vhat <- Matrix::solve(Vinv)
  
  Q_sp <- methods::as(methods::as(Q, "generalMatrix"), "CsparseMatrix")
  V_sp <- methods::as(methods::as(Vhat, "generalMatrix"), "CsparseMatrix")
  
  tr_cpp <- fast_trace_quad_cpp(Q_sp@i, Q_sp@p, Q_sp@x, V_sp@i, V_sp@p, V_sp@x)
  
  # Compare against S4 2-column index lookup and dense matrix trace
  idx_Q <- cbind(Q_sp@i + 1L, rep(1:ncol(Q_sp), diff(Q_sp@p)))
  tr_s4 <- sum(Q_sp@x * Vhat[idx_Q])
  tr_dense <- sum(diag(as.matrix(Q %*% Vhat)))
  
  expect_equal(tr_cpp, tr_s4, tolerance = 1e-12)
  expect_equal(tr_cpp, tr_dense, tolerance = 1e-12)
})
