# Form a sparse Nearest-Neighbor Gaussian Process Precision Matrix. See https://mc-stan.org/users/documentation/case-studies/nngp.html.

Form a sparse Nearest-Neighbor Gaussian Process Precision Matrix. See
https://mc-stan.org/users/documentation/case-studies/nngp.html.

## Usage

``` r
nngp_prec_mat(sp_dist, coords, cov_type, cov_params)
```

## Arguments

- sp_dist:

  Output of sparse_dist_LT function. Basically, top few rows are nearest
  distances, bottom few rows are indices. CALL THAT FUNCTION OR SEE IT
  FOR MORE DETAILS.

- coords:

  Matrix (x, y) of coordinates.

- cov_type:

  String for covariance model type. "Exp", "Sph", "Gau", and "Mat" are
  supported.

- cov_params:

  Covariance/kernel parameters. Same order as in variogram functions.
  Nugget variance (partial sill), Spatial variance (partial sill),
  Spatial range (scales distance); if Matern, smoothness kappa.

## Value

A list with a sparse precision matrix and associated eigenvalues.

Sparse precision matrix Q.

Eigenvalues Dinv.

Lower triangular factor A.

## Note

Applies to a single area.

Requires the Matrix library.

## Author

Florica J Constantine, florica AT berkeley.edu
