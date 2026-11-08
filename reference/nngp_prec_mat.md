# Form a sparse Nearest-Neighbor Gaussian Process Precision Matrix.

Form a sparse Nearest-Neighbor Gaussian Process Precision Matrix.

## Usage

``` r
nngp_prec_mat(sp_dist, nb_dist, cov_type, cov_params)
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

## Value

A list with a sparse precision matrix Q, eigenvalues Dinv, and factor A.

## Note

Applies to a single area.

## Author

Florica J Constantine, florica AT berkeley.edu
