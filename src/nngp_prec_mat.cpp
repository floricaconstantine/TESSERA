#include <Rcpp.h>
#include <Rmath.h>

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
//'
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
//' Performs the per-location covariance kernel evaluations and linear
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
//'
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
       
       for (int j = 0; j < n_keep; ++j) {
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
         if (i == j) {
           val_mat += nugget;
         }
         cov_nb_mat[i + j * n_keep] = val_mat;
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
     
     Dinv_val[idx] = 1.0 / (nugget + sill - dot_prod);
   }
   
   Dinv_val[n_cols] = 1.0 / (nugget + sill);
   
   return List::create(
     Named("A_row") = A_row,
     Named("A_col") = A_col,
     Named("A_val") = A_val,
     Named("Dinv")  = Dinv_val
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
//'
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
