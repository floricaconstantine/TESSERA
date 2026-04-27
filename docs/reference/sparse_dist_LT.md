# Get k nearest neighbors/distances given coordinates. Modified from: https://stackoverflow.com/questions/5560218/computing-sparse-pairwise-distance-matrix-in-r

Get k nearest neighbors/distances given coordinates. Modified from:
https://stackoverflow.com/questions/5560218/computing-sparse-pairwise-distance-matrix-in-r

## Usage

``` r
sparse_dist_LT(coords, k)
```

## Arguments

- coords:

  (x, y) for points as rows.

- k:

  Number of nearest neighbors.

## Value

2 k x n points matrix. First k rows are distances, last k are indices of
nearest neighbors. Only Lower Triangular part is formed.
