# M-step optimization for beta (covariate coefficients). Part of the M-Step in the EM algorithm.

M-step optimization for beta (covariate coefficients). Part of the
M-Step in the EM algorithm.

## Usage

``` r
M_step_beta(
  eta_list,
  Q_list,
  tau2_list,
  X_list,
  model_type = "spNNGP",
  gamma_list = NULL,
  XtDX_list = NULL,
  XtWX_list = NULL,
  XtWZX_list = NULL,
  XtX_list = NULL,
  XtDWIX_list = NULL
)
```

## Arguments

- eta_list:

  List of estimated means of eta.

- Q_list:

  List of unscaled precision matrices.

- tau2_list:

  List of precision/covariance matrix scalings.

- X_list:

  List of covariate matrices.

- model_type:

  "CAR", "SAR", "Leroux", or "spNNGP".

- gamma_list:

  List of current gamma estimates (required for lattice).

- XtDX_list:

  Precomputed X^T D X (required for CAR/SAR).

- XtWX_list:

  Precomputed X^T W X (required for CAR/SAR).

- XtWZX_list:

  Precomputed X^T W Z X (required for SAR).

- XtX_list:

  Precomputed X^T X (required for Leroux).

- XtDWIX_list:

  Precomputed X^T (D - W - I) X (required for Leroux).

## Value

Estimated beta vector.

## Author

Florica J Constantine, florica AT berkeley.edu
