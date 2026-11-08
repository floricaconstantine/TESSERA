# Maximize the expected likelihood in beta, holding other variables constant. Part of the M-Step in the EM algorithm.

Maximize the expected likelihood in beta, holding other variables
constant. Part of the M-Step in the EM algorithm.

## Usage

``` r
M_step_beta(
  eta_list,
  Q_list,
  tau2_list,
  X_list,
  model_type = "Leroux",
  gamma_list = NULL,
  XtDX_list = NULL,
  XtWX_list = NULL,
  XtWZX_list = NULL,
  XtX_list = NULL,
  XtDWIX_list = NULL,
  XtQX_list = NULL,
  XtQeta_list = NULL
)
```

## Arguments

- eta_list:

  List of estimated means of eta across all areas.

- Q_list:

  List of unscaled precision matrices across all areas.

- tau2_list:

  Vector of precision matrix scaling parameters tau^2 across all areas.

- X_list:

  List of covariate matrices across all areas.

- model_type:

  Character string specifying "CAR", "SAR", "Leroux", or "spNNGP".

- gamma_list:

  Vector of current spatial correlation parameters gamma across all
  areas.

- XtDX_list:

  Optional list of precomputed X^T D X matrices (CAR/SAR).

- XtWX_list:

  Optional list of precomputed X^T W X matrices (CAR/SAR).

- XtWZX_list:

  Optional list of precomputed X^T W Z X matrices (SAR).

- XtX_list:

  Optional list of precomputed X^T X matrices (Leroux).

- XtDWIX_list:

  Optional list of precomputed X^T (D - W - I) X matrices (Leroux).

- XtQX_list:

  Optional list of precomputed X^T Q X matrices (spNNGP inner CM steps).

- XtQeta_list:

  Optional list of precomputed X^T Q eta vectors (spNNGP inner CM
  steps).

## Value

Maximum likelihood estimate of covariate effects beta.

## Note

Applies across all areas.

Requires the Matrix library.

Requires the pracma library.

## Author

Florica J Constantine, florica AT berkeley.edu
