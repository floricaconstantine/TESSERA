# Fast C++ Evaluation of Lattice Precision Matrix Q from Pre-aligned CSC Vectors.

Computes Q = sum_m coeffs\[m\] \* mats\[m\] in a single O(nnz) loop and
returns a dgCMatrix sharing @i, @p, @Dim, and @Dimnames with the
precomputed template.

## Usage

``` r
eval_Q_lattice_cpp(precomp, coeffs)
```

## Arguments

- precomp:

  List returned by prep_Q_lattice_cpp.

- coeffs:

  NumericVector of scalar coefficients (length 2 or 3).

## Value

Updated precision matrix Q as a dgCMatrix.

## Author

Florica J Constantine, florica AT berkeley.edu
