# Raw Synthetic Spatial Data

A list containing raw components for three simulated spatial areas,
including counts, coordinates, and adjacency matrices.

## Format

A `.rds` file containing a list with 5 elements:

- count_matrix:

  A 1x900 matrix of simulated counts.

- meta_data:

  Data frame with sample identifiers.

- design_mat:

  Design matrix for fixed effects.

- coords:

  Matrix of spatial coordinates (x, y).

- W:

  Block-diagonal sparse adjacency matrix.
