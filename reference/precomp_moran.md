# Precompute coordinate-dependent spatial constants for Moran's I.

Since W, S1, S2, and EI depend solely on the spatial coordinates (c1,
c2) and not on the residual vector x, precomputing them once per area
avoids redundant work across EM iterations while preserving exact
output.

## Usage

``` r
precomp_moran(c1, c2)
```

## Arguments

- c1:

  x-coordinates.

- c2:

  y-coordinates.

## Value

NumericVector of length 4: c(W, S1, S2, ei).

## Author

Florica J Constantine, florica AT berkeley.edu
