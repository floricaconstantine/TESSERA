# Prepare data for the TESSERA method: Choose a distance threshold for adjacency matrices.

Prepare data for the TESSERA method: Choose a distance threshold for
adjacency matrices.

## Usage

``` r
visualizeNeighborDistances(meta_data, sample_col, coord_data, k_search = 20)
```

## Arguments

- meta_data:

  Dataframe with metadata/covariates.

- sample_col:

  String: Column name in meta_data identifying which rows correspond to
  which sample.

- coord_data:

  Dataframe or matrix with coordinates for observations.

- k_search:

  When forming adjacency matrices, the maximum number of neighbors.

## Value

A list comprised of distances and a ggplot2 visualization.

## Author

Florica J Constantine, florica AT berkeley.edu
