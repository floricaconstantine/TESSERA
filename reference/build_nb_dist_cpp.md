# Precompute Static Euclidean Distance Matrices for Neighbors in C++.

Replaces the R loop of as.matrix(stats::dist(...)) calls in
TESSERA_spNNGP.

## Usage

``` r
build_nb_dist_cpp(coords, sp_dist, nngp_k)
```

## Arguments

- coords:

  NumericMatrix of (x, y) coordinates (N x 2).

- sp_dist:

  NumericMatrix output of sparse_dist (2\*nngp_k x (N - 1)).

- nngp_k:

  Integer number of nearest neighbors.

## Value

A List of length ncol(sp_dist), where each element is the symmetric
Euclidean distance matrix among the valid neighbors of spot idx.

## Note

Applies to a single area.

## Author

Florica J Constantine, florica AT berkeley.edu
