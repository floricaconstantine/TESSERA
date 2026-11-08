# Compute unscaled precision matrix for SAR model.

Compute unscaled precision matrix for SAR model.

## Usage

``` r
Q_matrix_SAR(W, D, gamma_hat, precomp = NULL)
```

## Arguments

- W:

  Neighbor/adjacency matrix (symmetric, binary).

- D:

  Degree matrix (diagonal, values are row-sums of W).

- gamma_hat:

  Spatial correlation parameter.

- precomp:

  Optional precomputed list from `prep_Q_lattice_cpp` or precomputed
  matrix `W %*% D^{-1} %*% W`.

## Value

Unscaled precision matrix Q.

## Note

Applies to a single area.

Requires the Matrix library.

## Author

Florica J Constantine, florica AT berkeley.edu
