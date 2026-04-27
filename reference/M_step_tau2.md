# Maximize the expected likelihood in tau^2, holding other variables constant. Part of the M-Step in the EM algorithm.

Maximize the expected likelihood in tau^2, holding other variables
constant. Part of the M-Step in the EM algorithm.

## Usage

``` r
M_step_tau2(Vhat, eta_hat, Q, beta_hat, X)
```

## Arguments

- Vhat:

  Estimated covariance matrix of eta.

- eta_hat:

  Estimated mean of eta.

- Q:

  Unscaled precision matrix.

- beta_hat:

  Current estimate of covariate effects.

- X:

  Covariate matrix.

## Value

Maximum likelihood estimate of scaling parameter tau^2.

## Note

Applies to a single area.

## Author

Florica J Constantine, florica AT berkeley.edu
