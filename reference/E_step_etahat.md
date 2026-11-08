# Compute the Expectation of the random effects eta. Part of the E-Step in the EM algorithm. Follows Clayton and Kaldor (1987) for the Poisson-Gaussian model.

Compute the Expectation of the random effects eta. Part of the E-Step in
the EM algorithm. Follows Clayton and Kaldor (1987) for the
Poisson-Gaussian model.

## Usage

``` r
E_step_etahat(
  Vinv,
  Q,
  tau2,
  beta_hat,
  X,
  z,
  N,
  term1_precomp = NULL,
  chol_perm = NULL
)
```

## Arguments

- Vinv:

  Precomputed Precision matrix.

- Q:

  Unscaled precision matrix.

- tau2:

  Precision/covariance matrix scaling.

- beta_hat:

  Current estimate of covariate effects.

- X:

  Covariate matrix.

- z:

  Observed counts.

- N:

  Library size or expected counts.

- term1_precomp:

  Optional precomputed vector (z + 0.5)\*log((z + 0.5)/N) - 0.5.

- chol_perm:

  Optional list containing `Qpermchol` (lower-triangular Cholesky
  factor L) and `P` (permutation matrix) from
  `sparseinv::cholPermute(Vinv)` to solve via triangular substitution
  without re-factorizing `Vinv`.

## Value

Estimated mean of eta.

## Note

Applies to a single area.

This method is an Empirical Bayes approximation.

## Author

Florica J Constantine, florica AT berkeley.edu
