# Compute Moran's I using an adjacency matrix.

Compute Moran's I using an adjacency matrix.

## Usage

``` r
moran_I_nb(y, W)
```

## Arguments

- y:

  A vector of measurements.

- W:

  Neighbor/adjacency matrix (symmetric, binary).

## Value

Moran's I (spatial autocorrelation).

## Author

Florica J Constantine, florica AT berkeley.edu

## Examples

``` r
set.seed(2026)
tau2_true <- 1.0
gamma_true <- 0.5
beta_true <- c(1, 0, -1)
ex_data <- generate_data_one_area(1000, 0.03, "Leroux", beta_true,
  gamma_true, tau2_true, "rand_bern")
moran_I_nb(ex_data$z, ex_data$W)
#> [1] 0.1221835
```
