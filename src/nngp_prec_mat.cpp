#include <Rcpp.h>
#include <Rmath.h>
#include <algorithm>
#include <vector>
#include <string>
#include <limits>

using namespace Rcpp;

// Matern Kernel Function (C++ internal helper).
inline double eval_kernel_matern(double d, double sigma2, double rho, double nu, double const_term) {
  if (d == 0.0) {
    return sigma2;
  }
  double z = (std::sqrt(2.0 * nu) / rho) * std::abs(d);
  double z_term = nu * std::log(z);
  // R::bessel_k(x, nu, 1.0) matches base R besselK(x, nu, expon.scaled = FALSE)
  return std::exp(const_term + z_term) * R::bessel_k(z, nu, 1.0);
}

// Exponential Kernel Function (C++ internal helper).
inline double eval_kernel_exp(double d, double sigma2, double rho) {
  return sigma2 * std::exp(-std::abs(d / rho));
}

// Gaussian (Squared Exponential) Kernel Function (C++ internal helper).
inline double eval_kernel_gauss(double d, double sigma2, double rho) {
  double ratio = d / rho;
  return sigma2 * std::exp(-0.5 * (ratio * ratio));
}

// Spherical Kernel Function (C++ internal helper).
inline double eval_kernel_sph(double d, double sigma2, double rho) {
  double z = std::abs(d / rho);
  if (z > 1.0) {
    z = 1.0;
  }
  return sigma2 * (1.0 - 1.5 * z + 0.5 * (z * z * z));
}

// Candidate neighbor distance + 1-based index pair with deterministic tie-breaking
// matching base R's stable sort.list() behavior when distances are equal.
struct DistIdxPair {
  double dist;
  int idx;
  bool operator<(const DistIdxPair& other) const {
    if (dist != other.dist) {
      return dist < other.dist;
    }
    return idx < other.idx;
  }
};

//' Fast C++ Upper-Triangular Nearest Neighbor Distance Matrix.
//'
//' Computes the k nearest neighbors and Euclidean distances for each point i
//' among points with higher index (i + 1) .. n, matching sparse_dist in R.
//'
//' @author Florica J Constantine, florica AT berkeley.edu
//'
//' @note Applies to a single area.
//'
//' @param coords NumericMatrix of (x, y) coordinates (n x 2).
//' @param k Integer number of nearest neighbors.
//'
//' @returns NumericMatrix of dimension (2*k) x (n - 1), where rows 1..k contain
//'   Euclidean distances and rows (k + 1)..(2*k) contain 1-based neighbor indices.
// [[Rcpp::export]]
NumericMatrix sparse_dist_cpp(NumericMatrix coords, int k) {
 int n = coords.nrow();
 if (n < 2) {
   stop("sparse_dist_cpp: coords must have at least 2 rows.");
 }
 if (coords.ncol() < 2) {
   stop("sparse_dist_cpp: coords must have at least 2 columns.");
 }
 if (k < 1) {
   stop("sparse_dist_cpp: k must be >= 1.");
 }
 
 NumericMatrix out(2 * k, n - 1);
 CharacterVector rnames(2 * k);
 for (int r = 0; r < k; ++r) {
   rnames[r] = "d" + std::to_string(r + 1);
   rnames[k + r] = "i" + std::to_string(r + 1);
 }
 rownames(out) = rnames;
 
 std::vector<DistIdxPair> candidates;
 candidates.reserve(n);
 
 for (int i = 0; i < n - 1; ++i) {
   double xi = coords(i, 0);
   double yi = coords(i, 1);
   int n_cand = (n - 1) - i;
   candidates.resize(n_cand);
   
   for (int c = 0; c < n_cand; ++c) {
     int pt_idx = i + 1 + c;
     double dx = coords(pt_idx, 0) - xi;
     double dy = coords(pt_idx, 1) - yi;
     candidates[c].dist = std::sqrt(dx * dx + dy * dy);
     candidates[c].idx  = pt_idx + 1; // 1-based R index
   }
   
   int n_take = (k < n_cand) ? k : n_cand;
   std::partial_sort(candidates.begin(), candidates.begin() + n_take, candidates.end());
   
   for (int r = 0; r < n_take; ++r) {
     out(r, i)     = candidates[r].dist;
     out(k + r, i) = static_cast<double>(candidates[r].idx);
   }
   for (int r = n_take; r < k; ++r) {
     out(r, i)     = NA_REAL;
     out(k + r, i) = NA_REAL;
   }
 }
 
 return out;
}
 
//' Fast C++ Lower-Triangular Nearest Neighbor Distance Matrix.
//'
//' Computes the k nearest neighbors and Euclidean distances for each point i
//' among points with lower index 1 .. (i - 1), matching sparse_dist_LT in R.
//'
//' @author Florica J Constantine, florica AT berkeley.edu
//'
//' @note Applies to a single area.
//'
//' @param coords NumericMatrix of (x, y) coordinates (n x 2).
//' @param k Integer number of nearest neighbors.
//'
//' @returns NumericMatrix of dimension (2*k) x (n - 1), where rows 1..k contain
//'   Euclidean distances and rows (k + 1)..(2*k) contain 1-based neighbor indices.
// [[Rcpp::export]]
NumericMatrix sparse_dist_LT_cpp(NumericMatrix coords, int k) {
 int n = coords.nrow();
 if (n < 2) {
   stop("sparse_dist_LT_cpp: coords must have at least 2 rows.");
 }
 if (coords.ncol() < 2) {
   stop("sparse_dist_LT_cpp: coords must have at least 2 columns.");
 }
 if (k < 1) {
   stop("sparse_dist_LT_cpp: k must be >= 1.");
 }
 
 NumericMatrix out(2 * k, n - 1);
 CharacterVector rnames(2 * k);
 for (int r = 0; r < k; ++r) {
   rnames[r] = "d" + std::to_string(r + 1);
   rnames[k + r] = "i" + std::to_string(r + 1);
 }
 rownames(out) = rnames;
 
 std::vector<DistIdxPair> candidates;
 candidates.reserve(n);
 
 for (int i = 1; i < n; ++i) {
   double xi = coords(i, 0);
   double yi = coords(i, 1);
   int n_cand = i;
   candidates.resize(n_cand);
   
   for (int j = 0; j < n_cand; ++j) {
     double dx = coords(j, 0) - xi;
     double dy = coords(j, 1) - yi;
     candidates[j].dist = std::sqrt(dx * dx + dy * dy);
     candidates[j].idx  = j + 1; // 1-based R index
   }
   
   int n_take = (k < n_cand) ? k : n_cand;
   std::partial_sort(candidates.begin(), candidates.begin() + n_take, candidates.end());
   
   int col_out = i - 1;
   for (int r = 0; r < n_take; ++r) {
     out(r, col_out)     = candidates[r].dist;
     out(k + r, col_out) = static_cast<double>(candidates[r].idx);
   }
   for (int r = n_take; r < k; ++r) {
     out(r, col_out)     = NA_REAL;
     out(k + r, col_out) = NA_REAL;
   }
 }
 
 return out;
}
 
//' Precompute Static Euclidean Distance Matrices for Neighbors in C++.
//'
//' Replaces the R loop of as.matrix(stats::dist(...)) calls in TESSERA_spNNGP.
//'
//' @author Florica J Constantine, florica AT berkeley.edu
//'
//' @note Applies to a single area.
//'
//' @param coords NumericMatrix of (x, y) coordinates (N x 2).
//' @param sp_dist NumericMatrix output of sparse_dist (2*nngp_k x (N - 1)).
//' @param nngp_k Integer number of nearest neighbors.
//'
//' @returns A List of length ncol(sp_dist), where each element is the
//'   symmetric Euclidean distance matrix among the valid neighbors of spot idx.
// [[Rcpp::export]]
List build_nb_dist_cpp(NumericMatrix coords, NumericMatrix sp_dist, int nngp_k) {
 int n_cols = sp_dist.ncol();
 int n_rows_sp = sp_dist.nrow();
 int n_pts = coords.nrow();
 
 if (n_rows_sp < 2 * nngp_k) {
   stop("build_nb_dist_cpp: sp_dist has fewer than 2 * nngp_k rows.");
 }
 if (coords.ncol() < 2) {
   stop("build_nb_dist_cpp: coords must have at least 2 columns.");
 }
 
 List area_nb_dists(n_cols);
 
 for (int idx = 0; idx < n_cols; ++idx) {
   std::vector<int> keep_rows;
   keep_rows.reserve(nngp_k);
   for (int r = 0; r < nngp_k; ++r) {
     double nb_id = sp_dist(nngp_k + r, idx);
     if (!NumericVector::is_na(nb_id) && !ISNAN(nb_id)) {
       keep_rows.push_back(r);
     }
   }
   
   int n_keep = static_cast<int>(keep_rows.size());
   if (n_keep > 1) {
     NumericMatrix dmat(n_keep, n_keep);
     for (int i = 0; i < n_keep; ++i) {
       int row_i = static_cast<int>(sp_dist(nngp_k + keep_rows[i], idx)) - 1;
       if (row_i < 0 || row_i >= n_pts) {
         stop("build_nb_dist_cpp: neighbor index out of bounds.");
       }
       double xi = coords(row_i, 0);
       double yi = coords(row_i, 1);
       dmat(i, i) = 0.0;
       
       for (int j = i + 1; j < n_keep; ++j) {
         int row_j = static_cast<int>(sp_dist(nngp_k + keep_rows[j], idx)) - 1;
         if (row_j < 0 || row_j >= n_pts) {
           stop("build_nb_dist_cpp: neighbor index out of bounds.");
         }
         double dx = xi - coords(row_j, 0);
         double dy = yi - coords(row_j, 1);
         double d_ij = std::sqrt(dx * dx + dy * dy);
         dmat(i, j) = d_ij;
         dmat(j, i) = d_ij;
       }
     }
     area_nb_dists[idx] = dmat;
   } else {
     NumericMatrix dmat(1, 1);
     dmat(0, 0) = 0.0;
     area_nb_dists[idx] = dmat;
   }
 }
 
 return area_nb_dists;
}
 
//' Compute Sparse NNGP Precision Factors (A triplets and Dinv) in C++.
//'
//' Performs the per-location covariance kernel evaluations (exploiting symmetry
//' of the neighbor distance matrix to halve Bessel/exponential calls) and linear
//' system solves (LU decomposition with partial pivoting, matching base::solve)
//' for nngp_prec_mat.
//'
//' @author Florica J Constantine, florica AT berkeley.edu
//'
//' @note Applies to a single area.
//'
//' @param sp_dist NumericMatrix output of sparse_dist.
//' @param nb_dist List of precomputed Euclidean distance matrices for neighbors.
//' @param cov_type String for covariance model type ("Exp", "Sph", "Gau", "Mat").
//' @param cov_params NumericVector of covariance parameters in gstat convention:
//'   c(nugget, sill, range, \[smoothness\]).
//'
//' @returns A List containing:
//' \itemize{
//'   \item \strong{A_row}: 1-based row indices for sparse matrix A.
//'   \item \strong{A_col}: 1-based column indices for sparse matrix A.
//'   \item \strong{A_val}: Non-zero values for sparse matrix A.
//'   \item \strong{Dinv}: Vector of diagonal precision values D^\{-1\}.
//' }
// [[Rcpp::export]]
List nngp_prec_triplets_cpp(NumericMatrix sp_dist,
                           List nb_dist,
                           std::string cov_type,
                           NumericVector cov_params) {
 int nngp_k = sp_dist.nrow() / 2;
 int n_cols = sp_dist.ncol();
 
 if (cov_params.size() < 3) {
   stop("nngp_prec_triplets_cpp: cov_params must have length >= 3.");
 }
 if (nb_dist.size() < n_cols) {
   stop("nngp_prec_triplets_cpp: nb_dist length is smaller than ncol(sp_dist).");
 }
 
 double nugget = cov_params[0];
 double sill   = cov_params[1];
 double range  = cov_params[2];
 double nu     = (cov_params.size() >= 4 && !NumericVector::is_na(cov_params[3])) ? cov_params[3] : 0.5;
 
 int cov_code = 0;
 double matern_const = 0.0;
 if (cov_type == "Mat") {
   if (cov_params.size() < 4 || NumericVector::is_na(cov_params[3])) {
     stop("nngp_prec_triplets_cpp: Matern covariance requires a 4th parameter (nu/kappa).");
   }
   cov_code = 1;
   matern_const = std::log(sill) + (1.0 - nu) * std::log(2.0) - R::lgammafn(nu);
 } else if (cov_type == "Exp") {
   cov_code = 2;
 } else if (cov_type == "Gau") {
   cov_code = 3;
 } else if (cov_type == "Sph") {
   cov_code = 4;
 } else {
   stop("Invalid or unimplemented covariance structure.");
 }
 
 int total_nnz = 0;
 for (int idx = 0; idx < n_cols; ++idx) {
   for (int r = 0; r < nngp_k; ++r) {
     double nb_id = sp_dist(nngp_k + r, idx);
     if (!NumericVector::is_na(nb_id) && !ISNAN(nb_id)) {
       total_nnz++;
     }
   }
 }
 
 IntegerVector A_row(total_nnz);
 NumericVector A_col(total_nnz);
 NumericVector A_val(total_nnz);
 NumericVector Dinv_val(n_cols + 1, 0.0);
 
 std::vector<int> keep_idx;
 keep_idx.reserve(nngp_k);
 std::vector<double> cov_nb(nngp_k);
 std::vector<double> rhs(nngp_k);
 std::vector<double> cov_nb_mat(nngp_k * nngp_k);
 
 int pos = 0;
 double diag_val = sill + nugget;
 
 for (int idx = 0; idx < n_cols; ++idx) {
   keep_idx.clear();
   for (int r = 0; r < nngp_k; ++r) {
     double nb_id = sp_dist(nngp_k + r, idx);
     if (!NumericVector::is_na(nb_id) && !ISNAN(nb_id)) {
       keep_idx.push_back(r);
     }
   }
   
   int n_keep = static_cast<int>(keep_idx.size());
   if (n_keep == 0) {
     continue;
   }
   
   NumericMatrix nb_mat_in = as<NumericMatrix>(nb_dist[idx]);
   if (nb_mat_in.nrow() < n_keep || nb_mat_in.ncol() < n_keep) {
     stop("nngp_prec_triplets_cpp: nb_dist matrix dimensions smaller than n_keep.");
   }
   
   for (int i = 0; i < n_keep; ++i) {
     double d_vec = sp_dist(keep_idx[i], idx);
     double val_vec = 0.0;
     if (cov_code == 1) {
       val_vec = eval_kernel_matern(d_vec, sill, range, nu, matern_const);
     } else if (cov_code == 2) {
       val_vec = eval_kernel_exp(d_vec, sill, range);
     } else if (cov_code == 3) {
       val_vec = eval_kernel_gauss(d_vec, sill, range);
     } else {
       val_vec = eval_kernel_sph(d_vec, sill, range);
     }
     cov_nb[i] = val_vec;
     rhs[i]    = val_vec;
     
     // Diagonal entry has distance 0, so covariance is always sill + nugget
     cov_nb_mat[i + i * n_keep] = diag_val;
     
     // Evaluate only strict upper triangle j > i and mirror to (j, i)
     for (int j = i + 1; j < n_keep; ++j) {
       double d_ij = nb_mat_in(i, j);
       double val_mat = 0.0;
       if (cov_code == 1) {
         val_mat = eval_kernel_matern(d_ij, sill, range, nu, matern_const);
       } else if (cov_code == 2) {
         val_mat = eval_kernel_exp(d_ij, sill, range);
       } else if (cov_code == 3) {
         val_mat = eval_kernel_gauss(d_ij, sill, range);
       } else {
         val_mat = eval_kernel_sph(d_ij, sill, range);
       }
       cov_nb_mat[i + j * n_keep] = val_mat;
       cov_nb_mat[j + i * n_keep] = val_mat;
     }
   }
   
   // Solve cov_nb_mat * rhs = cov_nb via LU decomposition with partial pivoting
   for (int k = 0; k < n_keep; ++k) {
     int p = k;
     double max_val = std::abs(cov_nb_mat[k + k * n_keep]);
     for (int i = k + 1; i < n_keep; ++i) {
       double val = std::abs(cov_nb_mat[i + k * n_keep]);
       if (val > max_val) {
         max_val = val;
         p = i;
       }
     }
     if (max_val == 0.0 || ISNAN(max_val)) {
       stop("Singular neighbor covariance matrix in nngp_prec_triplets_cpp.");
     }
     if (p != k) {
       for (int j = 0; j < n_keep; ++j) {
         std::swap(cov_nb_mat[k + j * n_keep], cov_nb_mat[p + j * n_keep]);
       }
       std::swap(rhs[k], rhs[p]);
     }
     double pivot = cov_nb_mat[k + k * n_keep];
     for (int i = k + 1; i < n_keep; ++i) {
       double factor = cov_nb_mat[i + k * n_keep] / pivot;
       cov_nb_mat[i + k * n_keep] = factor;
       for (int j = k + 1; j < n_keep; ++j) {
         cov_nb_mat[i + j * n_keep] -= factor * cov_nb_mat[k + j * n_keep];
       }
     }
   }
   // Forward substitution (L * y = P * b, unit lower triangular L)
   for (int k = 0; k < n_keep; ++k) {
     for (int i = k + 1; i < n_keep; ++i) {
       rhs[i] -= cov_nb_mat[i + k * n_keep] * rhs[k];
     }
   }
   // Back substitution (U * x = y)
   for (int k = n_keep - 1; k >= 0; --k) {
     rhs[k] /= cov_nb_mat[k + k * n_keep];
     for (int i = 0; i < k; ++i) {
       rhs[i] -= cov_nb_mat[i + k * n_keep] * rhs[k];
     }
   }
   
   double dot_prod = 0.0;
   for (int i = 0; i < n_keep; ++i) {
     dot_prod += cov_nb[i] * rhs[i];
     A_row[pos] = idx + 2;
     A_col[pos] = sp_dist(nngp_k + keep_idx[i], idx);
     A_val[pos] = rhs[i];
     pos++;
   }
   
   Dinv_val[idx] = 1.0 / (diag_val - dot_prod);
 }
 
 Dinv_val[n_cols] = 1.0 / diag_val;
 
 return List::create(
   Named("A_row") = A_row,
   Named("A_col") = A_col,
   Named("A_val") = A_val,
   Named("Dinv")  = Dinv_val
 );
}
 
//' Compute X^T Q eta and X^T Q X Directly from NNGP Factors (A, Dinv) Without Forming Q.
//'
//' Since Q = (I - A)^T D^\{-1\} (I - A), we compute tilde_X = (I - A) X and
//' tilde_eta = (I - A) eta in O(nnz(A) * p) time, and then form
//' XtQeta = tilde_X^T D^\{-1\} tilde_eta and XtQX = tilde_X^T D^\{-1\} tilde_X.
//' This avoids building the sparse matrix Q during inner CM steps (opt_idx < opt_iters).
//'
//' @author Florica J Constantine, florica AT berkeley.edu
//'
//' @note Applies to a single area.
//'
//' @param A_row IntegerVector of 1-based row indices of A.
//' @param A_col NumericVector of 1-based column indices of A.
//' @param A_val NumericVector of non-zero values of A.
//' @param Dinv NumericVector of diagonal values of D^\{-1\} (length N).
//' @param X NumericMatrix of covariates (N x p).
//' @param eta NumericVector of random effect expectations (length N).
//'
//' @returns A List with:
//' \itemize{
//'   \item \strong{XtQeta}: NumericVector of length p equal to X^T Q eta.
//'   \item \strong{XtQX}: NumericMatrix of dimension (p x p) equal to X^T Q X.
//' }
// [[Rcpp::export]]
List nngp_beta_terms_cpp(const IntegerVector& A_row,
                        const NumericVector& A_col,
                        const NumericVector& A_val,
                        const NumericVector& Dinv,
                        const NumericMatrix& X,
                        const NumericVector& eta) {
 int N = Dinv.size();
 int p = X.ncol();
 int nnz = A_val.size();
 
 if (X.nrow() != N || eta.size() != N) {
   stop("nngp_beta_terms_cpp: X and eta dimensions must match length(Dinv).");
 }
 
 // Initialize tilde_eta = eta and tilde_X = X
 std::vector<double> tilde_eta(N);
 for (int i = 0; i < N; ++i) {
   tilde_eta[i] = eta[i];
 }
 
 std::vector<double> tilde_X(N * p);
 for (int m = 0; m < p; ++m) {
   int col_offset = m * N;
   for (int i = 0; i < N; ++i) {
     tilde_X[col_offset + i] = X(i, m);
   }
 }
 
 // Subtract A * eta and A * X in a single pass over the triplets of A
 for (int k = 0; k < nnz; ++k) {
   int r = A_row[k] - 1;
   int c = static_cast<int>(A_col[k]) - 1;
   double w = A_val[k];
   
   tilde_eta[r] -= w * eta[c];
   for (int m = 0; m < p; ++m) {
     tilde_X[m * N + r] -= w * X(c, m);
   }
 }
 
 // Compute XtQeta = tilde_X^T * diag(Dinv) * tilde_eta
 // and XtQX = tilde_X^T * diag(Dinv) * tilde_X
 NumericVector XtQeta(p);
 NumericMatrix XtQX(p, p);
 
 for (int m = 0; m < p; ++m) {
   int offset_m = m * N;
   double sum_eta = 0.0;
   for (int i = 0; i < N; ++i) {
     sum_eta += tilde_X[offset_m + i] * Dinv[i] * tilde_eta[i];
   }
   XtQeta[m] = sum_eta;
   
   for (int l = m; l < p; ++l) {
     int offset_l = l * N;
     double sum_X = 0.0;
     for (int i = 0; i < N; ++i) {
       sum_X += tilde_X[offset_m + i] * Dinv[i] * tilde_X[offset_l + i];
     }
     XtQX(m, l) = sum_X;
     XtQX(l, m) = sum_X;
   }
 }
 
 return List::create(
   Named("XtQeta") = XtQeta,
   Named("XtQX")   = XtQX
 );
}
 
//' Fast Trace of Product of Two Sparse CSC Matrices Tr(Q * Vhat).
//'
//' Computes sum(Q_sparse@x * Vhat\[idx_Q\], na.rm = TRUE) in O(nnz(Q) + nnz(V))
//' time without allocating a 2-column index matrix or invoking S4 subsetting.
//'
//' @author Florica J Constantine, florica AT berkeley.edu
//'
//' @param Q_i 0-based row indices (@i slot) of dgCMatrix Q.
//' @param Q_p Column pointers (@p slot) of dgCMatrix Q.
//' @param Q_x Non-zero values (@x slot) of dgCMatrix Q.
//' @param V_i 0-based row indices (@i slot) of dgCMatrix Vhat.
//' @param V_p Column pointers (@p slot) of dgCMatrix Vhat.
//' @param V_x Non-zero values (@x slot) of dgCMatrix Vhat.
//'
//' @returns Scalar trace value sum_\{i,j\} Q_\{i,j\} V_\{i,j\}.
// [[Rcpp::export]]
double fast_trace_quad_cpp(const IntegerVector& Q_i,
                          const IntegerVector& Q_p,
                          const NumericVector& Q_x,
                          const IntegerVector& V_i,
                          const IntegerVector& V_p,
                          const NumericVector& V_x) {
 int n_cols = Q_p.size() - 1;
 if (V_p.size() - 1 != n_cols) {
   stop("fast_trace_quad_cpp: Q and Vhat must have the same number of columns.");
 }
 
 double total = 0.0;
 
 for (int j = 0; j < n_cols; ++j) {
   int pQ = Q_p[j];
   int endQ = Q_p[j + 1];
   int pV = V_p[j];
   int endV = V_p[j + 1];
   
   while (pQ < endQ && pV < endV) {
     int rQ = Q_i[pQ];
     int rV = V_i[pV];
     if (rQ == rV) {
       double prod = Q_x[pQ] * V_x[pV];
       if (!ISNAN(prod)) {
         total += prod;
       }
       pQ++;
       pV++;
     } else if (rQ < rV) {
       pQ++;
     } else {
       pV++;
     }
   }
 }
 
 return total;
}
 
//' Pre-align Static Sparse Matrices onto a Shared CSC Template for Lattice Q Updates.
//'
//' Merges the sparsity patterns of 2 or 3 static dgCMatrix objects (such as
//' (D, W) for CAR, (I, D - W) for Leroux, or (D, W, W D^\{-1\} W) for SAR),
//' ensuring all diagonal entries are structurally present, and returns a
//' template dgCMatrix along with aligned numeric @x vectors for O(nnz) updates.
//'
//' @author Florica J Constantine, florica AT berkeley.edu
//'
//' @param mats List of 2 or 3 dgCMatrix objects of identical dimension (N x N).
//'
//' @returns A List with:
//' \itemize{
//'   \item \strong{tmpl}: Template dgCMatrix with merged @i and @p slots.
//'   \item \strong{x_list}: List of NumericVectors aligned to tmpl@i and tmpl@p.
//' }
// [[Rcpp::export]]
List prep_Q_lattice_cpp(List mats) {
 int K = mats.size();
 if (K < 2 || K > 3) {
   stop("prep_Q_lattice_cpp: mats must contain 2 or 3 dgCMatrix objects.");
 }
 
 S4 first_mat = as<S4>(mats[0]);
 IntegerVector dim_vec = first_mat.slot("Dim");
 int n_cols = dim_vec[1];
 
 std::vector<IntegerVector> i_vecs(K);
 std::vector<IntegerVector> p_vecs(K);
 std::vector<NumericVector> x_vecs(K);
 
 for (int m = 0; m < K; ++m) {
   S4 mat_m = as<S4>(mats[m]);
   i_vecs[m] = mat_m.slot("i");
   p_vecs[m] = mat_m.slot("p");
   x_vecs[m] = mat_m.slot("x");
 }
 
 std::vector<int> u_i;
 std::vector<int> u_p(n_cols + 1, 0);
 std::vector<std::vector<double>> aligned_x(K);
 
 std::vector<int> ptr(K);
 std::vector<int> end_ptr(K);
 
 for (int j = 0; j < n_cols; ++j) {
   for (int m = 0; m < K; ++m) {
     ptr[m] = p_vecs[m][j];
     end_ptr[m] = p_vecs[m][j + 1];
   }
   bool diag_done = false;
   
   while (true) {
     int min_r = std::numeric_limits<int>::max();
     for (int m = 0; m < K; ++m) {
       if (ptr[m] < end_ptr[m] && i_vecs[m][ptr[m]] < min_r) {
         min_r = i_vecs[m][ptr[m]];
       }
     }
     if (!diag_done && j < min_r) {
       min_r = j;
     }
     if (min_r == std::numeric_limits<int>::max()) {
       break;
     }
     if (min_r == j) {
       diag_done = true;
     }
     
     u_i.push_back(min_r);
     for (int m = 0; m < K; ++m) {
       if (ptr[m] < end_ptr[m] && i_vecs[m][ptr[m]] == min_r) {
         aligned_x[m].push_back(x_vecs[m][ptr[m]]);
         ptr[m]++;
       } else {
         aligned_x[m].push_back(0.0);
       }
     }
   }
   u_p[j + 1] = static_cast<int>(u_i.size());
 }
 
 int total_nnz = static_cast<int>(u_i.size());
 S4 tmpl("dgCMatrix");
 tmpl.slot("i") = wrap(u_i);
 tmpl.slot("p") = wrap(u_p);
 tmpl.slot("Dim") = dim_vec;
 tmpl.slot("Dimnames") = first_mat.slot("Dimnames");
 tmpl.slot("x") = NumericVector(total_nnz, 0.0);
 tmpl.slot("factors") = List();
 
 List x_list(K);
 for (int m = 0; m < K; ++m) {
   x_list[m] = wrap(aligned_x[m]);
 }
 
 return List::create(
   Named("tmpl")   = tmpl,
   Named("x_list") = x_list
 );
}
 
//' Fast C++ Evaluation of Lattice Precision Matrix Q from Pre-aligned CSC Vectors.
//'
//' Computes Q = sum_m coeffs\[m\] * mats\[m\] in a single O(nnz) loop and returns
//' a dgCMatrix sharing @i, @p, @Dim, and @Dimnames with the precomputed template.
//'
//' @author Florica J Constantine, florica AT berkeley.edu
//'
//' @param precomp List returned by prep_Q_lattice_cpp.
//' @param coeffs NumericVector of scalar coefficients (length 2 or 3).
//'
//' @returns Updated precision matrix Q as a dgCMatrix.
// [[Rcpp::export]]
S4 eval_Q_lattice_cpp(List precomp, NumericVector coeffs) {
 S4 tmpl = as<S4>(precomp["tmpl"]);
 List x_list = as<List>(precomp["x_list"]);
 int K = coeffs.size();
 
 NumericVector x0 = as<NumericVector>(x_list[0]);
 int nnz = x0.size();
 NumericVector out_x(nnz);
 
 if (K == 2) {
   double c0 = coeffs[0];
   double c1 = coeffs[1];
   NumericVector x1 = as<NumericVector>(x_list[1]);
   for (int idx = 0; idx < nnz; ++idx) {
     out_x[idx] = c0 * x0[idx] + c1 * x1[idx];
   }
 } else if (K == 3) {
   double c0 = coeffs[0];
   double c1 = coeffs[1];
   double c2 = coeffs[2];
   NumericVector x1 = as<NumericVector>(x_list[1]);
   NumericVector x2 = as<NumericVector>(x_list[2]);
   for (int idx = 0; idx < nnz; ++idx) {
     out_x[idx] = c0 * x0[idx] + c1 * x1[idx] + c2 * x2[idx];
   }
 } else {
   stop("eval_Q_lattice_cpp: coeffs must have length 2 or 3.");
 }
 
 S4 out("dgCMatrix");
 out.slot("i") = tmpl.slot("i");
 out.slot("p") = tmpl.slot("p");
 out.slot("Dim") = tmpl.slot("Dim");
 out.slot("Dimnames") = tmpl.slot("Dimnames");
 out.slot("x") = out_x;
 out.slot("factors") = List();
 
 return out;
}
 
//' Fast C++ Construction of E-Step Precision Matrix Vinv = Q / tau2 + diag(diag_add).
//'
//' Scales the non-zero entries of dgCMatrix Q by inv_tau2 and adds diag_add\[j\]
//' to each diagonal entry (j, j) without S4 Diagonal() coercions or sparse '+'.
//'
//' @author Florica J Constantine, florica AT berkeley.edu
//'
//' @param Q_dgC S4 dgCMatrix precision matrix Q (with all diagonal entries present).
//' @param inv_tau2 Scalar 1.0 / tau2.
//' @param diag_add NumericVector of diagonal increments (0.5 + z).
//'
//' @returns S4 dgCMatrix Vinv sharing @i, @p, @Dim, and @Dimnames with Q_dgC.
// [[Rcpp::export]]
S4 build_Vinv_cpp(S4 Q_dgC, double inv_tau2, NumericVector diag_add) {
 IntegerVector Q_i = Q_dgC.slot("i");
 IntegerVector Q_p = Q_dgC.slot("p");
 NumericVector Q_x = Q_dgC.slot("x");
 
 int n_cols = Q_p.size() - 1;
 int nnz = Q_x.size();
 NumericVector V_x(nnz);
 
 for (int j = 0; j < n_cols; ++j) {
   double d_add = diag_add[j];
   int start = Q_p[j];
   int end   = Q_p[j + 1];
   for (int p = start; p < end; ++p) {
     double val = Q_x[p] * inv_tau2;
     if (Q_i[p] == j) {
       val += d_add;
     }
     V_x[p] = val;
   }
 }
 
 S4 Vinv("dgCMatrix");
 Vinv.slot("i") = Q_i;
 Vinv.slot("p") = Q_p;
 Vinv.slot("Dim") = Q_dgC.slot("Dim");
 Vinv.slot("Dimnames") = Q_dgC.slot("Dimnames");
 Vinv.slot("x") = V_x;
 Vinv.slot("factors") = List();
 
 return Vinv;
}
 