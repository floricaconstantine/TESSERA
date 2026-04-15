# Fast computation of Moran's I.

Fast computation of Moran's I.

## Usage

``` r
calc_moran(x, c1, c2)
```

## Arguments

- x:

  Measurements.

- c1:

  x-coordinates.

- c2:

  y-coordinates.

## Value

I Moran's I.

EI Expected value of Moran's I under the null of no spatial correlation.

SD Standard deviation of Moran's I.

## Author

Matthew Cooper.

## Examples

``` r
calc_moran(rpois(1000, 1), runif(1000), runif(1000))
#> [1] -0.001938802 -0.001001001  0.002866237
```
