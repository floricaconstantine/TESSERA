# Compute the unscaled precision matrix in a Leroux model.

Compute the unscaled precision matrix in a Leroux model.

## Usage

``` r
Q_matrix_Leroux(W, D, gamma_val)
```

## Arguments

- W:

  Neighbor/adjacency matrix (symmetric, binary).

- D:

  Degree matrix (diagonal, values are row-sums of W).

- gamma_val:

  Correlation parameter.

## Value

Unscaled precision matrix.

## Note

Applies to a single area.

Requires the Matrix library.

## Author

Florica J Constantine, florica AT berkeley.edu
