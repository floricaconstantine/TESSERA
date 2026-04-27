# Compute the log likelihood of a set of Poisson variables.

Compute the log likelihood of a set of Poisson variables.

## Usage

``` r
poisson_loglike(z, theta_hat)
```

## Arguments

- z:

  Observed counts.

- theta_hat:

  Estimated Poisson parameters.

## Value

log likelihood of data.

## Note

Discards invalid values in computation.

Discards factorial terms for numerical stability, etc.

## Author

Florica J Constantine, florica AT berkeley.edu
