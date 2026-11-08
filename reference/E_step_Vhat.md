# Compute the Covariance Matrix of the random effects eta. Part of the E-Step in the EM algorithm. Follows Clayton and Kaldor (1987) for the Poisson-Gaussian model.

Compute the Covariance Matrix of the random effects eta. Part of the
E-Step in the EM algorithm. Follows Clayton and Kaldor (1987) for the
Poisson-Gaussian model.

## Usage

``` r
E_step_Vhat(Vinv, P = NULL, cholQp = NULL, return_dgC = FALSE)
```

## Arguments

- Vinv:

  Precomputed Precision matrix.

- P:

  Cholesky Permutation (precomputed).

- cholQp:

  Optional precomputed lower-triangular Cholesky factor of
  `t(P) %*% Vinv %*% P` from
  [`sparseinv::cholPermute`](https://rdrr.io/pkg/sparseinv/man/cholPermute.html).

- return_dgC:

  Logical: If `TRUE`, returns the sparse inverse subset directly as a
  `dgCMatrix` without coercing to `symmetricMatrix`. Defaults to `FALSE`
  for backward compatibility.

## Value

Estimated covariance matrix (sparse subset).

## Note

Applies to a single area.

This method is an Empirical Bayes approximation.

This method computes a sparse inverse subset using the Takahashi
equations to avoid an O(N^3) dense matrix inversion. It only computes
the diagonal and non-zero off-diagonal elements that match the sparsity
pattern of the precision matrix.

Requires the Matrix library.

Requires the sparseinv library.

## Author

Florica J Constantine, florica AT berkeley.edu
