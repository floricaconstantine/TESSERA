# Compute the expected log likelihood.

Compute the expected log likelihood.

## Usage

``` r
expected_loglike(
  Vhat,
  eta_hat,
  Q,
  gamma_hat,
  tau2,
  beta_hat,
  X,
  W,
  D,
  eig_vals,
  model_type
)
```

## Arguments

- Vhat:

  Estimated covariance matrix of eta.

- eta_hat:

  Estimated mean of eta.

- Q:

  Unscaled precision matrix. The scaled precision matrix in spNNGP.

- gamma_hat:

  Parameter of interest: Estimate of correlation parameter.

- tau2:

  Parameter of interest: Current estimate of scale.

- beta_hat:

  Current estimate of covariate effects.

- X:

  Covariate matrix.

- W:

  Neighbor/adjacency matrix (symmetric, binary).

- D:

  Degree matrix (diagonal, values are row-sums of W).

- eig_vals:

  Eigenvalues of Z = D^{-1} W (CAR/SAR), D - W (Leroux), or Q (spNNGP).

- model_type:

  "CAR", "SAR", or "Leroux", or "spNNGP". If spNNGP, set tau^2 to 1 and
  the value of gamma is irrelevant. Model for the random effects.

## Value

Expected log likelihood.

## Note

Applies to a single area.

Multivariate normal; drops the 2 pi term.

## Author

Florica J Constantine, florica AT berkeley.edu
