# Maximize the expected likelihood in gamma, holding other variables constant. Part of the M-Step in the EM algorithm. ONLY APPLIES TO THE SAR MODEL.

Maximize the expected likelihood in gamma, holding other variables
constant. Part of the M-Step in the EM algorithm. ONLY APPLIES TO THE
SAR MODEL.

## Usage

``` r
M_step_gamma_SAR(
  trace_scalars,
  eta_hat,
  tau2,
  beta_hat,
  X,
  W,
  D,
  eig_vals,
  gamma_current = NULL
)
```

## Arguments

- trace_scalars:

  List of precomputed trace scalars from the E-step.

- eta_hat:

  Estimated mean of eta.

- tau2:

  Precision matrix scaling parameter tau^2.

- beta_hat:

  Current estimate of covariate effects.

- X:

  Covariate matrix.

- W:

  Neighbor/adjacency matrix (symmetric, binary).

- D:

  Diagonal degree matrix (row sums of W).

- eig_vals:

  Eigenvalues of Z = D^{-1} W.

- gamma_current:

  Previous value of gamma. Only needed for failure modes.

## Value

Maximum likelihood estimate of correlation parameter gamma. A list with
gamma_hat (estimate) and grad_val (gradient).

## Note

Applies to a single area.

Requires the Matrix library.

Requires the pracma library.

## Author

Florica J Constantine, florica AT berkeley.edu
