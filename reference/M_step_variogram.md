# Holding beta constant, fit a variogram model to eta - X beta. Part of the M-Step in the EM algorithm. Instead of a traditional MLE, which has cubic time complexity, we fit a variogram to estimate the kernel parameters, which is quadratic time. ONLY APPLIES TO THE GAUSSIAN PROCESS MODEL.

Holding beta constant, fit a variogram model to eta - X beta. Part of
the M-Step in the EM algorithm. Instead of a traditional MLE, which has
cubic time complexity, we fit a variogram to estimate the kernel
parameters, which is quadratic time. ONLY APPLIES TO THE GAUSSIAN
PROCESS MODEL.

## Usage

``` r
M_step_variogram(eta_hat, beta_hat, X, coords, cov_type = "Exp")
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

## Value

A vector of fitted parameters. Nugget variance (partial sill), Spatial
variance (partial sill), Spatial range (scales distance); if Matern,
smoothness kappa.

## Note

Applies to a single area.

Requires the sp library.

Requires the gstat library.

## Author

Florica J Constantine, florica AT berkeley.edu
