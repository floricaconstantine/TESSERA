# Compute the Negative Hessian of gamma. ONLY APPLIES TO THE CAR MODEL.

Compute the Negative Hessian of gamma. ONLY APPLIES TO THE CAR MODEL.

## Usage

``` r
neg_hessian_gamma_CAR(gamma_hat, eig_vals)
```

## Arguments

- gamma_hat:

  Parameter of interest: Estimate of correlation parameter.

- eig_vals:

  Eigenvalues of Z = D^{-1} W.

## Value

Negative Hessian of gamma. Also called the observed information.

## Note

Applies to a single area.

## Author

Florica J Constantine, florica AT berkeley.edu
