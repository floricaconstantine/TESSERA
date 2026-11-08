# Fast C++ Lower-Triangular Nearest Neighbor Distance Matrix.

Computes the k nearest neighbors and Euclidean distances for each point
i among points with lower index 1 .. (i - 1), matching sparse_dist_LT in
R.

## Usage

``` r
sparse_dist_LT_cpp(coords, k)
```

## Arguments

- coords:

  NumericMatrix of (x, y) coordinates (n x 2).

- k:

  Integer number of nearest neighbors.

## Value

NumericMatrix of dimension (2*k) x (n - 1), where rows 1..k contain
Euclidean distances and rows (k + 1)..(2*k) contain 1-based neighbor
indices.

## Note

Applies to a single area.

## Author

Florica J Constantine, florica AT berkeley.edu
