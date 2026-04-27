# Optimal threshold selection for empirical null estimation

In generalized linear mixed models (GLMMs) and the TESSERA framework,
theoretical null distributions for Wald statistics are often unreliable
due to finite-sample biases and E-step approximations. This function
estimates an empirical null distribution—modeled as a scaled,
non-central \\\chi^2_1\\ distribution—by identifying an optimal
threshold below which statistics likely originate from the null.

## Usage

``` r
select_Wald_threshold(wald_stats, quantile_spacing = 0.01, metric = "Raw_MSE")
```

## Arguments

- wald_stats:

  A numeric vector of Wald statistics (non-negative). Note: t-statistics
  should be squared to obtain F-statistics (or 1-df \\\chi^2\\
  equivalent) before input.

- quantile_spacing:

  Numeric value defining the search resolution for the threshold.
  Defaults to 0.01 (searching 0.01, 0.02, ..., 0.99 quantiles).

- metric:

  Character string specifying the selection metric. Options include
  "Raw\_" or "Log\_" prefixed to "MSE", "MAE", "MedAE", or "MaxAE".
  These correspond to Mean Squared Error, Mean Absolute Error, Median
  Absolute Error, and Maximum Absolute Error (L-infinity), respectively.

## Value

A list containing the following components:

- **threshold**: The selected optimal Wald statistic threshold.

- **chi2params**: A vector containing the estimated shift
  (non-centrality) and scale parameters for the empirical null.

- **threshold_results**: A data frame summarizing metrics and fitted
  parameters across all evaluated threshold quantiles.

## Details

The estimation procedure sweeps over a range of candidate thresholds.
For each threshold, it estimates parameters such that the resulting
p-values best approximate a \\Uniform(0, 1)\\ distribution. The optimal
threshold is selected by minimizing a chosen error metric (e.g., MSE)
between the observed p-value quantiles and theoretical uniform
quantiles.

## Author

Florica J Constantine, florica AT berkeley.edu

## Examples

``` r
# Simulate null and non-null Wald statistics
null_stats <- 5 * stats::rchisq(2000, df = 1, ncp = 0.5)
alt_stats <- 200 + 5 * stats::rchisq(500, df = 1, ncp = 10)
wald_vec <- c(null_stats, alt_stats)

# Select optimal threshold
select_Wald_threshold(wald_vec, quantile_spacing = 0.05)
#> Error in `$<-`(`*tmp*`, "quantile", value = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)): Assigned data `quantile_list` must be compatible with existing data.
#> ✖ Existing data has 0 rows.
#> ✖ Assigned data has 19 rows.
#> ℹ Only vectors of size 1 are recycled.
#> Caused by error in `vectbl_recycle_rhs_rows()`:
#> ! Can't recycle input of size 19 to size 0.
```
