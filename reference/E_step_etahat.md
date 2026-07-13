# Compute the Expectation of the random effects eta. Part of the E-Step in the EM algorithm. Follows Clayton and Kaldor (1987) for the Poisson-Gaussian model.

Compute the Expectation of the random effects eta. Part of the E-Step in
the EM algorithm. Follows Clayton and Kaldor (1987) for the
Poisson-Gaussian model.

## Usage

``` r
E_step_etahat(Q, tau2, beta_hat, X, z, N)
```

## Arguments

- Q:

  Unscaled precision matrix.

- tau2:

  Precision/covariance matrix scaling.

- beta_hat:

  Current estimate of covariate effects.

- X:

  Covariate matrix.

- z:

  Observed counts.

- N:

  Library size or expected counts.

## Value

Estimated mean of eta.

## Note

Applies to a single area.

This method is an Empirical Bayes approximation.

## Author

Florica J Constantine, florica AT berkeley.edu
