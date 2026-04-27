# Compute the Negative Hessian of gamma. ONLY APPLIES TO THE SAR MODEL.

Compute the Negative Hessian of gamma. ONLY APPLIES TO THE SAR MODEL.

## Usage

``` r
neg_hessian_gamma_SAR(
  Vhat,
  eta_hat,
  gamma_hat,
  tau2,
  beta_hat,
  X,
  W,
  D,
  eig_vals
)
```

## Arguments

- Vhat:

  Estimated covariance matrix of eta.

- eta_hat:

  Estimated mean of eta.

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

  Eigenvalues of Z = D^{-1} W.

## Value

Negative Hessian of gamma. Also called the observed information.

## Note

Applies to a single area.

Requires the Matrix library.

## Author

Florica J Constantine, florica AT berkeley.edu
