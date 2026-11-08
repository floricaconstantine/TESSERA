# Compute X^T Q eta and X^T Q X Directly from NNGP Factors (A, Dinv) Without Forming Q.

Since Q = (I - A)^T D^{-1} (I - A), we compute tilde_X = (I - A) X and
tilde_eta = (I - A) eta in O(nnz(A) \* p) time, and then form XtQeta =
tilde_X^T D^{-1} tilde_eta and XtQX = tilde_X^T D^{-1} tilde_X. This
avoids building the sparse matrix Q during inner CM steps (opt_idx \<
opt_iters).

## Usage

``` r
nngp_beta_terms_cpp(A_row, A_col, A_val, Dinv, X, eta)
```

## Arguments

- A_row:

  IntegerVector of 1-based row indices of A.

- A_col:

  NumericVector of 1-based column indices of A.

- A_val:

  NumericVector of non-zero values of A.

- Dinv:

  NumericVector of diagonal values of D^{-1} (length N).

- X:

  NumericMatrix of covariates (N x p).

- eta:

  NumericVector of random effect expectations (length N).

## Value

A List with:

- **XtQeta**: NumericVector of length p equal to X^T Q eta.

- **XtQX**: NumericMatrix of dimension (p x p) equal to X^T Q X.

## Note

Applies to a single area.

## Author

Florica J Constantine, florica AT berkeley.edu
