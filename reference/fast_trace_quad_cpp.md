# Fast Trace of Product of Two Sparse CSC Matrices Tr(Q \* Vhat).

Computes sum(Q_sparse@x \* Vhat\[idx_Q\], na.rm = TRUE) in O(nnz(Q) +
nnz(V)) time without allocating a 2-column index matrix or invoking S4
subsetting.

## Usage

``` r
fast_trace_quad_cpp(Q_i, Q_p, Q_x, V_i, V_p, V_x)
```

## Arguments

- Q_i:

  0-based row indices (@i slot) of dgCMatrix Q.

- Q_p:

  Column pointers (@p slot) of dgCMatrix Q.

- Q_x:

  Non-zero values (@x slot) of dgCMatrix Q.

- V_i:

  0-based row indices (@i slot) of dgCMatrix Vhat.

- V_p:

  Column pointers (@p slot) of dgCMatrix Vhat.

- V_x:

  Non-zero values (@x slot) of dgCMatrix Vhat.

## Value

Scalar trace value sum\_{i,j} Q\_{i,j} V\_{i,j}.

## Author

Florica J Constantine, florica AT berkeley.edu
