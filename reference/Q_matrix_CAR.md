# Compute unscaled precision matrix for CAR model.

Compute unscaled precision matrix for CAR model.

## Usage

``` r
Q_matrix_CAR(W, D, gamma_hat, precomp = NULL)
```

## Arguments

- W:

  Neighbor/adjacency matrix (symmetric, binary).

- D:

  Degree matrix (diagonal, values are row-sums of W).

- gamma_hat:

  Spatial correlation parameter.

- precomp:

  Optional precomputed list from `prep_Q_lattice_cpp`.

## Value

Unscaled precision matrix Q.

## Note

Applies to a single area.

## Author

Florica J Constantine, florica AT berkeley.edu
