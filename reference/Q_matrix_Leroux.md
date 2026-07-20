# Compute the unscaled precision matrix in a Leroux model.

Compute the unscaled precision matrix in a Leroux model.

## Usage

``` r
Q_matrix_Leroux(W, D, gamma_val, precomp = NULL)
```

## Arguments

- W:

  Neighbor/adjacency matrix (symmetric, binary).

- D:

  Degree matrix (diagonal, values are row-sums of W).

- gamma_val:

  Correlation parameter.

- precomp:

  Precomputed list containing D_minus_W and id_mat.

## Value

Unscaled precision matrix.

## Note

Applies to a single area.

Requires the Matrix library.

## Author

Florica J Constantine, florica AT berkeley.edu
