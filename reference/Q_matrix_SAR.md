# Compute the unscaled precision matrix in a SAR model.

Compute the unscaled precision matrix in a SAR model.

## Usage

``` r
Q_matrix_SAR(W, D, gamma_val)
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
