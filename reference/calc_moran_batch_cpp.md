# Batch Computation of Moran's I Across Multiple Vectors on Shared Coordinates.

Computes the Moran's I statistic (first element of calc_moran) for all M
columns of X_mat in a single upper-triangular pass over (c1, c2). Used
by summarize_TESSERA to compute all 8 summary Moran's I metrics
simultaneously.

## Usage

``` r
calc_moran_batch_cpp(X_mat, c1, c2)
```

## Arguments

- X_mat:

  NumericMatrix of dimension (N x M), where each column is a measurement
  vector on the same N spatial locations.

- c1:

  NumericVector of x-coordinates (length N).

- c2:

  NumericVector of y-coordinates (length N).

## Value

NumericVector of length M containing Moran's I for each column.

## Author

Florica J Constantine, florica AT berkeley.edu
