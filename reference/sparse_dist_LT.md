# Get k nearest neighbors/distances given coordinates.

For each point, we only look at points with a lower index. E.g., for
point 10, we only look at points 1, 2, ..., 9. We do not look at the
first point (since there are no lower-index points). So the output has
one fewer column than the number of points. If there are fewer than k
points with lower index, we fill the remaining entries with NA.

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

2 k x n points matrix. Top k rows: Euclidean distances to nearest
neighbors. Bottom k rows: Indices of nearest neighbors. Each column
corresponds to a point.

## Note

Applies to a single area.

## Author

Florica J Constantine, florica AT berkeley.edu

## Examples

``` r
sparse_dist_LT(matrix(1:4, ncol=2), 1)
#>        [,1]
#> d1 1.414214
#> i1 1.000000
```
