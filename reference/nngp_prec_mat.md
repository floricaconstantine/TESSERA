# Form a sparse Nearest-Neighbor Gaussian Process Precision Matrix.

Form a sparse Nearest-Neighbor Gaussian Process Precision Matrix.

## Usage

``` r
nngp_prec_mat(
  sp_dist,
  nb_dist,
  cov_type,
  cov_params,
  build_Q = TRUE,
  X = NULL,
  eta = NULL
)
```

## Arguments

- sp_dist:

  Output of sparse_dist_LT function.

- nb_dist:

  Precomputed list of Euclidean distance matrices for neighbors.

- cov_type:

  String for covariance model type.

- cov_params:

  Covariance/kernel parameters. Assuming gstat conventions, nugget
  (non-spatial), sill (spatial/scale), range, order (Matern).

- build_Q:

  Logical: If `TRUE` (default), constructs the sparse precision matrix
  `Q = (I - A)^T D^{-1} (I - A)` and sparse factor `A`. If `FALSE`,
  skips building `Q` and `A` and instead computes `XtQeta` and `XtQX`
  directly in C++ from `X` and `eta`.

- X:

  Optional covariate matrix (N x p) used when `build_Q = FALSE`.

- eta:

  Optional random effect vector (length N) used when `build_Q = FALSE`.

## Value

A list with a sparse precision matrix Q, eigenvalues Dinv, and factor A
(plus `XtQeta` and `XtQX` when `build_Q = FALSE`).

## Note

Applies to a single area.

## Author

Florica J Constantine, florica AT berkeley.edu
