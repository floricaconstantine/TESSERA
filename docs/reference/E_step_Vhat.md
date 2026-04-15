# Compute the Covariance Matrix of the random effects eta. Part of the E-Step in the EM algorithm. Follows Clayton and Kaldor (1987) for the Poisson-Gaussian model.

Compute the Covariance Matrix of the random effects eta. Part of the
E-Step in the EM algorithm. Follows Clayton and Kaldor (1987) for the
Poisson-Gaussian model.

## Usage

``` r
E_step_Vhat(Q, tau2, z)
```

## Arguments

- Q:

  Unscaled precision matrix.

- tau2:

  Precision/covariance matrix scaling.

- z:

  Observed counts.

## Value

Estimated covariance matrix.

## Note

Applies to a single area.

This method is an Empirical Bayes approximation.

This method inverts a matrix, because the inverse is actually needed...

Requires the Matrix library.

Requires the pracma library.

## Author

Florica J Constantine, florica AT berkeley.edu
