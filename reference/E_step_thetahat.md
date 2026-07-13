# Compute the estimated poisson parameters theta. Note that eta = log theta. Associated with the E-Step in the EM algorithm, but not needed to run. I.e., this is a utility function. Follows Clayton and Kaldor (1987) for the Poisson-Gaussian model.

Compute the estimated poisson parameters theta. Note that eta = log
theta. Associated with the E-Step in the EM algorithm, but not needed to
run. I.e., this is a utility function. Follows Clayton and Kaldor (1987)
for the Poisson-Gaussian model.

## Usage

``` r
E_step_thetahat(Vhat, eta_hat)
```

## Arguments

- Vhat:

  Estimated covariance matrix of eta (sparse subset).

- eta_hat:

  Estimated mean of eta.

## Value

Estimated poisson parameters theta_hat.

## Note

Applies to a single area.

## Author

Florica J Constantine, florica AT berkeley.edu
