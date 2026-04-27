# Compute the INVERSE Covariance Matrix of the random effects eta. Part of the E-Step in the EM algorithm. Follows Clayton and Kaldor (1987) for the Poisson-Gaussian model.

Compute the INVERSE Covariance Matrix of the random effects eta. Part of
the E-Step in the EM algorithm. Follows Clayton and Kaldor (1987) for
the Poisson-Gaussian model.

## Usage

``` r
E_step_Vhat_PLU(Q, tau2, z)
```

## Arguments

- Q:

  Unscaled precision matrix.

- tau2:

  Precision/covariance matrix scaling.

- z:

  Observed counts.

## Value

Estimated covariance matrix P LU decomposition.

## Note

Applies to a single area.

This method is an Empirical Bayes approximation.

Requires the Matrix library.

## Author

Florica J Constantine, florica AT berkeley.edu
