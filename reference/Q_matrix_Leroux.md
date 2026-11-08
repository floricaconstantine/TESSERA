# Compute unscaled precision matrix for Leroux model.

Compute unscaled precision matrix for Leroux model.

## Usage

``` r
Q_matrix_Leroux(W, D, gamma_hat, precomp = NULL)
```

## Arguments

- W:

  Neighbor/adjacency matrix (symmetric, binary).

- D:

  Degree matrix (diagonal, values are row-sums of W).

- gamma_hat:

  Spatial correlation parameter.

- precomp:

  Optional precomputed list from `prep_Q_lattice_cpp` or list containing
  `id_mat` and `D_minus_W`.

## Value

Unscaled precision matrix Q.

## Note

Applies to a single area.

Requires the Matrix library.

## Author

Florica J Constantine, florica AT berkeley.edu
