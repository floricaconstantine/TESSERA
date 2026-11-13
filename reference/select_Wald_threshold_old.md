# Optimal threshold selection for empirical null estimation

In generalized linear mixed models (GLMMs) and the TESSERA framework,
theoretical null distributions for Wald statistics are often unreliable
due to finite-sample biases and E-step approximations. This function
estimates an empirical null distribution—modeled as a scaled,
non-central \\\chi^2_1\\ distribution—by identifying an optimal
threshold below which statistics likely originate from the null.

## Usage

``` r
select_Wald_threshold_old(
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
select_Wald_threshold(wald_vec, quantile_spacing = 0.05)
#> $threshold
#> [1] 121.3377
#> 
#> $chi2_params
#>     scale     shift 
#> 4.7788778 0.6473782 
#> 
#> $threshold_results
#>      threshold       scale        shift n_obs      Raw_MSE     Raw_MAE
#> 1    0.0658458  0.01492689 1.423003e+00   125 0.0115709042 0.092201225
#> 2    0.2476369  1.02838733 0.000000e+00   250 0.1298542210 0.311308745
#> 3    0.4950019  1.03720314 2.000000e+00   375 0.2040149686 0.391477936
#> 4    0.8541667  0.37055072 1.592741e+00   500 0.0538145748 0.201063960
#> 5    1.3743795  0.41899265 1.545149e+00   625 0.0273402557 0.143220652
#> 6    2.0311153  1.04136254 1.213334e+00   750 0.0508386287 0.194621932
#> 7    2.9267893  0.99782664 1.210663e+00   875 0.0248341119 0.135878743
#> 8    3.8939586  3.46914457 7.300028e-01  1000 0.0664964937 0.222742241
#> 9    5.2985028  2.33224446 8.549648e-01  1125 0.0276035325 0.143114046
#> 10   6.8360407  5.88629421 3.434530e-01  1250 0.0431456160 0.179397407
#> 11   8.9845325  3.82580357 6.202439e-01  1375 0.0189507342 0.118655258
#> 12  11.2754844  9.41896446 1.559298e-03  1500 0.0251865733 0.137451740
#> 13  14.3048693  9.93382982 9.348559e-06  1625 0.0177646747 0.115691709
#> 14  17.7699534 11.27580972 0.000000e+00  1750 0.0148680267 0.106731413
#> 15  26.2685437  4.61888471 6.838310e-01  1875 0.0012270363 0.030790716
#> 16 121.3377326  4.77887782 6.473782e-01  2000 0.0000306778 0.004684714
#> 17 231.3650773 20.53941093 9.203559e-07  2125 0.0140173466 0.105009449
#> 18 250.7221020 34.62095347 1.387459e-06  2250 0.0268229801 0.143073043
#> 19 278.4050579 51.23475190 1.058365e-07  2375 0.0339819487 0.158318098
#>      Raw_MedAE  Raw_MaxAE      Log_MSE     Log_MAE  Log_MedAE Log_MaxAE
#> 1  0.084273324 0.18126565 0.0866300149 0.163968333 0.06764605 1.6657349
#> 2  0.300679671 0.62431106 0.2588118887 0.339019864 0.20449219 2.4951593
#> 3  0.394190853 0.78170689 0.3101272693 0.382899249 0.25246022 2.7688453
#> 4  0.192556249 0.40082600 0.1815871920 0.265736202 0.14124346 2.6040380
#> 5  0.142768425 0.28601344 0.1363609455 0.217624423 0.10891256 2.5539699
#> 6  0.196565829 0.39011024 0.1774338674 0.261094159 0.14399254 2.7680201
#> 7  0.144462060 0.27268719 0.1311462319 0.211156588 0.11005914 2.6796118
#> 8  0.228350301 0.44620761 0.1989760448 0.282092497 0.16337052 2.9510533
#> 9  0.148713482 0.28719139 0.1376142371 0.218052544 0.11308292 2.8110254
#> 10 0.183545762 0.35931222 0.1663276683 0.249086872 0.13576149 2.9538951
#> 11 0.117787666 0.23802787 0.1174703754 0.194558348 0.09102109 2.8166234
#> 12 0.142736534 0.27476561 0.1321690762 0.212318905 0.10906305 2.9153241
#> 13 0.117938718 0.23029879 0.1137240615 0.190843820 0.09140116 2.8744430
#> 14 0.110512644 0.21022490 0.1041076549 0.180478597 0.08672467 2.8649722
#> 15 0.032161913 0.06711708 0.0344821715 0.078496637 0.02783730 2.3544515
#> 16 0.004561014 0.01282857 0.0003479877 0.007965561 0.00443816 0.5491459
#> 17 0.116217511 0.17830686 0.1332229009 0.195711963 0.10422135 1.5576521
#> 18 0.155385086 0.25360502 0.0871472872 0.205672030 0.14082517 1.5058624
#> 19 0.165655265 0.29432007 0.0618838817 0.191228126 0.15575287 1.9722631
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
