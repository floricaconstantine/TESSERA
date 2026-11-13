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
#> $threshold
#> [1] 109.9827
#> 
#> $chi2_params
#>     scale     shift 
#> 4.5074751 0.7227001 
#> 
#> $threshold_results
#>       threshold       scale        shift n_obs      Raw_MSE     Raw_MAE
#> 1    0.04325045  0.01847568 1.197236e+00   125 3.773261e-02 0.167779801
#> 2    0.20142869  0.30574554 0.000000e+00   250 5.784954e-02 0.207502489
#> 3    0.51129626  1.06610854 3.869237e-02   375 8.281518e-02 0.251967147
#> 4    0.96692831  1.79505726 0.000000e+00   500 7.168426e-02 0.234088550
#> 5    1.48372394  2.00000000 7.353031e-01   625 9.768686e-02 0.272991747
#> 6    2.11613750  2.00000000 1.066886e+00   750 9.034651e-02 0.262883534
#> 7    2.99037463  2.62293234 8.621395e-01   875 7.289809e-02 0.235351006
#> 8    4.03990795  3.33138199 7.625583e-01  1000 6.281763e-02 0.218223438
#> 9    5.40075626  2.40907065 9.040613e-01  1125 2.993309e-02 0.151008291
#> 10   7.00458394  5.99118544 3.844595e-01  1250 4.486984e-02 0.184170751
#> 11   8.56689053  7.64216080 4.685971e-01  1375 5.202394e-02 0.199169928
#> 12  10.64388631  9.17332878 3.815748e-01  1500 4.525118e-02 0.185907365
#> 13  13.63181917  5.69317411 6.111268e-01  1625 1.793213e-02 0.117045473
#> 14  17.65773282  4.56965906 7.669544e-01  1750 6.544241e-03 0.071209190
#> 15  25.55272203  3.83425924 9.065845e-01  1875 9.063396e-04 0.027096393
#> 16 109.98273467  4.50747506 7.227001e-01  2000 3.157371e-05 0.004398976
#> 17 227.53104835 20.29810726 2.053179e-07  2125 1.389384e-02 0.104082958
#> 18 248.52970428 34.19703376 1.212065e-06  2250 2.672200e-02 0.142381423
#> 19 271.29183235 50.80154161 1.447859e-07  2375 3.407493e-02 0.158157497
#>     Raw_MedAE  Raw_MaxAE      Log_MSE     Log_MAE   Log_MedAE Log_MaxAE
#> 1  0.16262940 0.33796192 0.1520117154 0.238355895 0.119771455 1.9319177
#> 2  0.20824569 0.42171507 0.1872848941 0.270724616 0.151216347 2.3227687
#> 3  0.25396526 0.49562215 0.2144699411 0.301012133 0.178381330 2.5713787
#> 4  0.23873398 0.46246817 0.2042024633 0.289105443 0.169377160 2.6660199
#> 5  0.27097194 0.54053727 0.2321270739 0.315249186 0.188068569 2.8303643
#> 6  0.26291261 0.52075174 0.2247659330 0.308659841 0.183534633 2.8932776
#> 7  0.22956201 0.46701802 0.2069481917 0.290417954 0.163901781 2.9129027
#> 8  0.21015332 0.43367563 0.1947404844 0.278344341 0.152342473 2.9385643
#> 9  0.14651784 0.30144321 0.1428633723 0.224296935 0.111610509 2.8294372
#> 10 0.18109308 0.36648667 0.1691204737 0.252510399 0.133960705 2.9624718
#> 11 0.19968329 0.39443678 0.1792290027 0.263802288 0.145931500 3.0357102
#> 12 0.19083036 0.36770009 0.1691884970 0.253383829 0.140361369 3.0430084
#> 13 0.12022691 0.23534160 0.1153179395 0.191995571 0.093339291 2.8779898
#> 14 0.07361443 0.14252593 0.0739937812 0.139510812 0.059058930 2.6907148
#> 15 0.02548666 0.05581728 0.0297226277 0.070793159 0.023266421 2.2896394
#> 16 0.00355186 0.01621876 0.0003685056 0.009117609 0.003284448 0.3262052
#> 17 0.11236500 0.18344534 0.1341429703 0.195647309 0.100974577 1.5887179
#> 18 0.15148073 0.25448001 0.0874719230 0.205528894 0.138960591 1.4996305
#> 19 0.16235801 0.29533604 0.0621800371 0.191200792 0.153486552 1.9957306
#>    quantile
#> 1      0.05
#> 2      0.10
#> 3      0.15
#> 4      0.20
#> 5      0.25
#> 6      0.30
#> 7      0.35
#> 8      0.40
#> 9      0.45
#> 10     0.50
#> 11     0.55
#> 12     0.60
#> 13     0.65
#> 14     0.70
#> 15     0.75
#> 16     0.80
#> 17     0.85
#> 18     0.90
#> 19     0.95
#> 
```
