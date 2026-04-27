# Compute the Negative Hessian for beta.

Compute the Negative Hessian for beta.

## Usage

``` r
neg_hessian_beta(Q_list, tau2_list, X_list)
```

## Arguments

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

Negative Hessian matrix of beta. Also called the observed information
matrix.

## Note

Applies to ALL areas.

Requires the Matrix library.

## Author

Florica J Constantine, florica AT berkeley.edu
