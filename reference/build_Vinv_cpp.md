# Fast C++ Construction of E-Step Precision Matrix Vinv = Q / tau2 + diag(diag_add).

Scales the non-zero entries of dgCMatrix Q by inv_tau2 and adds
diag_add\[j\] to each diagonal entry (j, j) without S4 Diagonal()
coercions or sparse '+'.

## Usage

``` r
build_Vinv_cpp(Q_dgC, inv_tau2, diag_add)
```

## Arguments

- Q_dgC:

  S4 dgCMatrix precision matrix Q (with all diagonal entries present).

- inv_tau2:

  Scalar 1.0 / tau2.

- diag_add:

  NumericVector of diagonal increments (0.5 + z).

## Value

S4 dgCMatrix Vinv sharing @i, @p, @Dim, and @Dimnames with Q_dgC.

## Author

Florica J Constantine, florica AT berkeley.edu
