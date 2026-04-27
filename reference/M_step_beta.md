# Maximize the expected likelihood in beta, holding other variables constant. Part of the M-Step in the EM algorithm.

Maximize the expected likelihood in beta, holding other variables
constant. Part of the M-Step in the EM algorithm.

## Usage

``` r
M_step_beta(eta_list, Q_list, tau2_list, X_list)
```

## Arguments

- eta_list:

  List of the estimated means of eta. One vector per area.

- Q_list:

  List of the unscaled precision matrices. One matrix per area, same
  ordering and length as eta_list. Depends on gamma, so that dependency
  is implicit. I.e., Q_list should be updated in each iteration after
  gamma is updated.

- tau2_list:

  List or vector of precision matrix scaling values. Same length and
  ordering as eta_list, with one number per area.

- X_list:

  List of covariate matrices. One matrix per area, same ordering and
  length as eta_list.

## Value

Maximum likelihood estimate of covariates beta.

## Note

Applies to ALL areas.

Requires the Matrix library.

Requires the pracma library.

## Author

Florica J Constantine, florica AT berkeley.edu
