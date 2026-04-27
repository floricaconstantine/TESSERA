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
#> [1] 109.9827
#> 
#> $chi2_params
#>     scale     shift 
#> 4.5074748 0.7227003 
#> 
#> $threshold_results
#>        threshold       scale        shift n_obs      Raw_MSE     Raw_MAE
#> 5%    0.04325045  0.01847573 1.197234e+00   125 3.773271e-02 0.167780006
#> 10%   0.20142869  0.06258233 1.120973e+00   250 1.815473e-02 0.117183546
#> 15%   0.51129626  1.11207752 4.053210e-07   375 8.308833e-02 0.252378564
#> 20%   0.96692831  1.79505558 9.206677e-07   500 7.168426e-02 0.234088548
#> 25%   1.48372394  2.00000000 7.353034e-01   625 9.768688e-02 0.272991778
#> 30%   2.11613750  2.00000000 1.066886e+00   750 9.034650e-02 0.262883524
#> 35%   2.99037463  2.62293234 8.621391e-01   875 7.289806e-02 0.235350961
#> 40%   4.03990795  3.33137732 7.625596e-01  1000 6.281760e-02 0.218223400
#> 45%   5.40075626  2.40905735 9.040661e-01  1125 2.993297e-02 0.151008011
#> 50%   7.00458394  5.99118544 3.844596e-01  1250 4.486985e-02 0.184170764
#> 55%   8.56689053  7.64216080 4.685975e-01  1375 5.202397e-02 0.199169982
#> 60%  10.64388631  9.17332878 3.815745e-01  1500 4.525116e-02 0.185907332
#> 65%  13.63181917  5.69316743 6.111277e-01  1625 1.793211e-02 0.117045407
#> 70%  17.65773282  4.56965856 7.669547e-01  1750 6.544245e-03 0.071209211
#> 75%  25.55272203  3.83425812 9.065849e-01  1875 9.063388e-04 0.027096385
#> 80% 109.98273467  4.50747477 7.227003e-01  2000 3.157368e-05 0.004398972
#> 85% 227.53104835 20.29812611 0.000000e+00  2125 1.389386e-02 0.104083058
#> 90% 248.52970428 34.19724618 0.000000e+00  2250 2.672222e-02 0.142381970
#> 95% 271.29183235 50.80155825 0.000000e+00  2375 3.407494e-02 0.158157512
#>       Raw_MedAE  Raw_MaxAE      Log_MSE     Log_MAE   Log_MedAE Log_MaxAE
#> 5%  0.162629623 0.33796239 0.1520118966 0.238356075 0.119771566 1.9319183
#> 10% 0.119647434 0.24595470 0.1145899987 0.192077735 0.094291897 2.0746654
#> 15% 0.254387273 0.49640503 0.2147613564 0.301289604 0.178624349 2.5720623
#> 20% 0.238733981 0.46246817 0.2042024621 0.289105442 0.169377159 2.6660199
#> 25% 0.270971976 0.54053732 0.2321270942 0.315249205 0.188068588 2.8303644
#> 30% 0.262912597 0.52075173 0.2247659263 0.308659835 0.183534626 2.8932775
#> 35% 0.229561964 0.46701794 0.2069481619 0.290417923 0.163901751 2.9129026
#> 40% 0.210153292 0.43367553 0.1947404503 0.278344311 0.152342454 2.9385642
#> 45% 0.146517602 0.30144252 0.1428631039 0.224296670 0.111610349 2.8294362
#> 50% 0.181093092 0.36648669 0.1691204829 0.252510409 0.133960715 2.9624719
#> 55% 0.199683353 0.39443687 0.1792290402 0.263802329 0.145931537 3.0357103
#> 60% 0.190830320 0.36770003 0.1691884745 0.253383804 0.140361347 3.0430083
#> 65% 0.120226846 0.23534146 0.1153178797 0.191995502 0.093339247 2.8779895
#> 70% 0.073614456 0.14252594 0.0739937904 0.139510833 0.059058951 2.6907148
#> 75% 0.025486653 0.05581725 0.0297226144 0.070793137 0.023266418 2.2896392
#> 80% 0.003551874 0.01621873 0.0003685057 0.009117608 0.003284431 0.3262052
#> 85% 0.112365105 0.18344551 0.1341427343 0.195647305 0.100974662 1.5887162
#> 90% 0.151481265 0.25448109 0.0874710442 0.205528566 0.138961042 1.4996392
#> 95% 0.162358028 0.29533608 0.0621800248 0.191200781 0.153486566 1.9957308
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
