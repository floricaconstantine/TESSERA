# Given the output of the TESSERA algorithms and a contrast matrix, compute Wald T-statistics.

Given the output of the TESSERA algorithms and a contrast matrix,
compute Wald T-statistics.

## Usage

``` r
calc_Wald_statistics_old(TESSERAOutput_obj, contrast_mat)
```

## Arguments

- TESSERAOutput_obj:

  Output of the TESSERA algorithms.

- contrast_mat:

  Matrix with contrasts of the estimated coeficients. Rows are
  contrasts, columns correspond to columns in beta_hat (covariates).
  beta_hat is the vector of estimated coefficients, stored in
  TESSERAOutput_obj.

## Value

A dataframe of Wald t-statistics. Square these to get F-statistics.

## Author

Florica J Constantine, florica AT berkeley.edu

## Examples

``` r
# Locate the saved model results in inst/extdata
rds_path <- system.file("extdata", "example_TESSERA_out_Leroux.rds", package = "TESSERA")
# Load the results object
TESSERA_out_Leroux <- readRDS(rds_path)

# Run the Wald test
calc_Wald_statistics(TESSERAOutput_obj = TESSERA_out_Leroux,
                 contrast_mat = matrix(c(1, 0, 0), nrow = 1))
#>      gene fit_model kernel_type contrast_string contrast_indices contrast_val
#> 1 example    Leroux          NA        1*+0*+0*                1     1.021903
#>   contrast_se wald_stat_t
#> 1  0.05245812    19.48037
```
