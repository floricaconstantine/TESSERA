# Compute the best fit to the data. Associated with the E-Step in the EM algorithm, but not needed to run. I.e., this is a utility function. Follows Clayton and Kaldor (1987) for the Poisson-Gaussian model.

Compute the best fit to the data. Associated with the E-Step in the EM
algorithm, but not needed to run. I.e., this is a utility function.
Follows Clayton and Kaldor (1987) for the Poisson-Gaussian model.

## Usage

``` r
E_step_predict(theta_hat, N)
```

## Arguments

- theta_hat:

  Estimated poisson parameters.

- N:

  Library size or expected counts.

## Value

Estimated values z (model fits).

## Author

Florica J Constantine, florica AT berkeley.edu
