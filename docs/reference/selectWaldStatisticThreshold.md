# Optimal threshold selection for empirical null estimation

In generalized linear mixed models (GLMMs) and the TESSERA framework,
theoretical null distributions for Wald statistics are often unreliable
due to finite-sample biases and E-step approximations. This function
estimates an empirical null distribution—modeled as a scaled,
non-central \\\chi^2_1\\ distribution—by identifying an optimal
threshold below which statistics likely originate from the null.

## Usage

``` r
selectWaldStatisticThreshold(
  wald_stats,
  quantile_spacing = 0.01,
  metric = "Raw_MSE"
)
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
selectWaldStatisticThreshold(wald_vec, quantile_spacing = 0.05)
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> $threshold
#> [1] 103.8969
#> 
#> $chi2_params
#>     scale     shift 
#> 5.2776758 0.4643671 
#> 
#> $threshold_results
#>        threshold      scale        shift n_obs      Raw_MSE     Raw_MAE
#> 5%    0.04972148  0.1830294 3.890485e-07   125 1.213025e-01 0.302174062
#> 10%   0.19191037  0.1555879 1.218828e+00   250 8.702800e-02 0.256127124
#> 15%   0.49428686  0.1387604 1.203979e+00   375 1.560678e-02 0.108727612
#> 20%   0.90655979  1.8338922 1.013005e-05   500 7.750195e-02 0.241500833
#> 25%   1.46457814  0.9246736 7.707580e-01   625 4.518213e-02 0.184795387
#> 30%   2.09631064  2.0000000 7.852391e-01   750 7.463585e-02 0.238043855
#> 35%   2.78918973  2.5687734 1.073974e+00   875 8.866900e-02 0.259884428
#> 40%   3.75489422  3.3664869 9.347429e-01  1000 7.862041e-02 0.244319361
#> 45%   4.94593628  4.2808722 7.440705e-01  1125 6.522067e-02 0.221999957
#> 50%   6.32611289  3.5232420 8.484486e-01  1250 4.073378e-02 0.175616039
#> 55%   7.98650758  5.2114239 6.478470e-01  1375 4.155790e-02 0.177132267
#> 60%  10.65189858  2.7752883 9.663906e-01  1500 9.208840e-03 0.083879969
#> 65%  13.53805125  4.1316246 7.060520e-01  1625 9.668258e-03 0.085429976
#> 70%  17.70162169  5.3801723 5.015304e-01  1750 6.580981e-03 0.070160329
#> 75%  26.65311100  3.9824434 7.560453e-01  1875 6.335301e-04 0.022139291
#> 80% 103.89687505  5.2776758 4.643671e-01  2000 2.169594e-05 0.003681268
#> 85% 228.21603609 20.2396930 0.000000e+00  2125 1.494732e-02 0.107985158
#> 90% 246.60636185 34.1445590 0.000000e+00  2250 2.781201e-02 0.145348213
#> 95% 269.61658037 50.7555100 0.000000e+00  2375 3.502884e-02 0.160487225
#>      Raw_MedAE  Raw_MaxAE      Log_MSE     Log_MAE   Log_MedAE Log_MaxAE
#> 5%  0.30558068 0.59990260 0.2479219181 0.332671815 0.207375166 2.1789069
#> 10% 0.25303075 0.50976395 0.2202766074 0.304349323 0.177846039 2.4080397
#> 15% 0.10050899 0.21648187 0.1065180711 0.183357852 0.079311664 2.2102953
#> 20% 0.24655203 0.48144235 0.2113819314 0.294833058 0.174090913 2.6833605
#> 25% 0.18866009 0.37132024 0.1698925692 0.252917800 0.139034910 2.6637186
#> 30% 0.23483282 0.47348509 0.2082741926 0.292221720 0.167139609 2.8517592
#> 35% 0.26252795 0.51615156 0.2232275607 0.306852957 0.183306433 2.9562958
#> 40% 0.24174267 0.48556804 0.2129665751 0.296596223 0.171218618 2.9877270
#> 45% 0.21969692 0.44219834 0.1978905793 0.281200976 0.158179640 2.9982358
#> 50% 0.16898375 0.34962130 0.1628790369 0.245560203 0.126029218 2.9416960
#> 55% 0.17468162 0.35317180 0.1640403630 0.246898496 0.130128874 2.9871658
#> 60% 0.07800327 0.16744021 0.0854877760 0.155583359 0.062710147 2.6980919
#> 65% 0.08382924 0.17068935 0.0870997076 0.157731767 0.067315835 2.7425941
#> 70% 0.07031575 0.14015990 0.0731301865 0.139105507 0.057852591 2.6915761
#> 75% 0.02191422 0.04635234 0.0247409550 0.061784174 0.018571162 2.2098495
#> 80% 0.00314378 0.01522793 0.0002265485 0.006159724 0.003539694 0.4360268
#> 85% 0.11796584 0.18219258 0.1344090136 0.197439707 0.109205726 1.5783751
#> 90% 0.15450258 0.25418059 0.0876306857 0.206883176 0.144872012 1.5107005
#> 95% 0.16300190 0.29351720 0.0622222549 0.192103009 0.159016096 2.0025864
#>     quantile
#> 5%      0.05
#> 10%     0.10
#> 15%     0.15
#> 20%     0.20
#> 25%     0.25
#> 30%     0.30
#> 35%     0.35
#> 40%     0.40
#> 45%     0.45
#> 50%     0.50
#> 55%     0.55
#> 60%     0.60
#> 65%     0.65
#> 70%     0.70
#> 75%     0.75
#> 80%     0.80
#> 85%     0.85
#> 90%     0.90
#> 95%     0.95
#> 
```
