# Pre-align Static Sparse Matrices onto a Shared CSC Template for Lattice Q Updates.

Merges the sparsity patterns of 2 or 3 static dgCMatrix objects (such as
(D, W) for CAR, (I, D - W) for Leroux, or (D, W, W D^{-1} W) for SAR),
ensuring all diagonal entries are structurally present, and returns a
template dgCMatrix along with aligned numeric @x vectors for O(nnz)
updates.

## Usage

``` r
prep_Q_lattice_cpp(mats)
```

## Arguments

- mats:

  List of 2 or 3 dgCMatrix objects of identical dimension (N x N).

## Value

A List with:

- **tmpl**: Template dgCMatrix with merged @i and @p slots.

- **x_list**: List of NumericVectors aligned to tmpl@i and tmpl@p.

## Author

Florica J Constantine, florica AT berkeley.edu
