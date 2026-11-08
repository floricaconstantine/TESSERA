# Compute Sparse NNGP Precision Factors (A triplets and Dinv) in C++.

Performs the per-location covariance kernel evaluations and linear
system solves (LU decomposition with partial pivoting, matching
base::solve) for nngp_prec_mat.

## Usage

``` r
nngp_prec_triplets_cpp(sp_dist, nb_dist, cov_type, cov_params)
```

## Arguments

- sp_dist:

  NumericMatrix output of sparse_dist.

- nb_dist:

  List of precomputed Euclidean distance matrices for neighbors.

- cov_type:

  String for covariance model type ("Exp", "Sph", "Gau", "Mat").

- cov_params:

  NumericVector of covariance parameters in gstat convention: c(nugget,
  sill, range, \[smoothness\]).

## Value

A List containing:

- **A_row**: 1-based row indices for sparse matrix A.

- **A_col**: 1-based column indices for sparse matrix A.

- **A_val**: Non-zero values for sparse matrix A.

- **Dinv**: Vector of diagonal precision values D^{-1}.

## Note

Applies to a single area.

## Author

Florica J Constantine, florica AT berkeley.edu
