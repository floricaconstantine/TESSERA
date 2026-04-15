# Holding beta constant, use BRISC to find the spatial parameters of eta - X beta. Part of the M-Step in the EM algorithm. ONLY APPLIES TO THE GAUSSIAN PROCESS MODEL.

Holding beta constant, use BRISC to find the spatial parameters of eta -
X beta. Part of the M-Step in the EM algorithm. ONLY APPLIES TO THE
GAUSSIAN PROCESS MODEL.

## Usage

``` r
M_step_BRISC(eta_hat, beta_hat, X, coords, cov_type = "Exp", k = 15)
```

## Arguments

- eta_hat:

  Estimated mean of eta.

- beta_hat:

  Current estimate of covariate effects.

- X:

  Covariate matrix.

- coords:

  Matrix of (x, y) coordinates.

- cov_type:

  String for covariance model type. "Exp", "Sph", "Gau", and "Mat" are
  supported.

- k:

  Number of neighbors.

## Value

A vector of fitted parameters. Nugget variance (partial sill), Spatial
variance (partial sill), Spatial range (scales distance); if Matern,
smoothness kappa.

## Note

Applies to a single area.

Requires the BRISC library (called indirectly).

## Author

Florica J Constantine, florica AT berkeley.edu
