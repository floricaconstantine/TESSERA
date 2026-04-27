# Wrapper to fit sparse NN GP model via BRISC.

This function theoretically works for multiple areas by stacking
everything together and ignoring differences in coordinates.

## Usage

``` r
BRISC_wrapper(
  z_list,
  X_list,
  coords_list,
  k = 15,
  cov_type = "Exp",
  transform_z = TRUE,
  z_offset = 0.5,
  verbose = FALSE
)
```

## Arguments

- z_list:

  List of vectors of observed counts—one vector per area.

- X_list:

  List of design or covariate matrices—one matrix per area. Same length
  and ordering as z_list. Matrices with number of rows equal to length
  of corresponding vector in z_list.

- coords_list:

  List of coordinate matrices (x, y)—one matrix per area. Same length
  and ordering as z_list. Matrices with number of rows equal to length
  of corresponding vector in z_list.

- k:

  Number of neighbors.

- cov_type:

  String for covariance model type. "Exp", "Sph", "Gau", and "Mat" are
  supported.

- transform_z:

  Boolean: log-transform z or not.

- z_offset:

  If transform_z, the counts z are transformed as log(z + z_offset).

- verbose:

  Whether to print output as BRISC runs.

## Value

The output list from BRISC, plus a time field for how long the function
ran for.

## Note

Requires the BRISC library.

## Author

Florica J Constantine, florica AT berkeley.edu
