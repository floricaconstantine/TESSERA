# Fit a scaled non-central chi^2_1 distribution using BOBYQA

Numerically maximizes the likelihood of a scaled non-central chi-square
distribution (df=1) truncated above a specific threshold. The fit is
conditional on the threshold.

## Usage

``` r
fit_scaled_noncentral_chi2_old(wald_stats, wald_thresh)
```

## Arguments

- wald_stats:

  Vector of Wald statistics (non-negative).

- wald_thresh:

  Threshold below which data is included in the fit.

## Value

Vector of c(scale, shift).

## Examples

``` r
fit_scaled_noncentral_chi2(5 * stats::rchisq(1000, 1, ncp=10), 20)
#> scaling   shift 
#>      NA      NA 
```
