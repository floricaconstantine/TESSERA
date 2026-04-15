# Compute the Negative Hessian for tau^2.

Compute the Negative Hessian for tau^2.

## Usage

``` r
neg_hessian_tau2(Vhat, eta_hat, Q, tau2_hat, beta_hat, X)
```

## Arguments

- Vhat:

  Estimated covariance matrix of eta.

- eta_hat:

  Estimated mean of eta.

- Q:

  Unscaled precision matrix.

- tau2_hat:

  Parameter of interest: Current estimate of scale.

- beta_hat:

  Current estimate of covariate effects.

- X:

  Covariate matrix.

## Value

Negative Hessian of tau^2. Also called the observed information.

## Note

Applies to a single area.

## Author

Florica J Constantine, florica AT berkeley.edu
