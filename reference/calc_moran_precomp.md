# Fast computation of Moran's I using precomputed coordinate constants.

Returns the exact same c(I, ei, sd) vector as calc_moran, using the
precomputed W, S1, S2, and ei from precomp_moran and an upper-triangular
pass over pairs (i \< j).

## Usage

``` r
calc_moran_precomp(x, c1, c2, moran_consts)
```

## Arguments

- x:

  Measurements (e.g., residuals).

- c1:

  x-coordinates.

- c2:

  y-coordinates.

- moran_consts:

  NumericVector of length 4 from precomp_moran: c(W, S1, S2, ei).

## Value

NumericVector of length 3: c(I, ei, sd).

## Author

Florica J Constantine, florica AT berkeley.edu
