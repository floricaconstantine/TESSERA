# The TESSERAOutput Class

The `TESSERAOutput` object is the primary model-fitting output container
for the `TESSERA` package, returned by
[`TESSERA_lattice`](https://floricaconstantine.github.io/TESSERA/reference/TESSERA_lattice.md)
and
[`TESSERA_spNNGP`](https://floricaconstantine.github.io/TESSERA/reference/TESSERA_spNNGP.md).
It stores estimated fixed effects, sample-specific spatial random
effects, spatial covariance parameters, optimization histories,
goodness-of-fit diagnostics, and negative Hessian matrices for
downstream inference.

## Details

Downstream functions—-specifically
[`calc_Wald_statistics`](https://floricaconstantine.github.io/TESSERA/reference/calc_Wald_statistics.md)
(for linear contrast testing) and
[`summarize_TESSERA`](https://floricaconstantine.github.io/TESSERA/reference/summarize_TESSERA.md)
(for sample-level performance summaries)—-ingest `TESSERAOutput` objects
directly.

## Elements

A `TESSERAOutput` object is a structured list containing the following
components:

- **beta_hat**: Numeric vector of estimated shared fixed effect
  coefficients \\\hat{\boldsymbol{\beta}}\\.

- **gamma_hat**: Numeric vector of estimated spatial
  correlation/dependence parameters \\\hat{\gamma}\_i\\ for each sample
  (*lattice models only*).

- **tau2_hat**: Numeric vector of estimated spatial variance scale
  parameters \\\hat{\tau}\_i^2\\ for each sample (*lattice models
  only*).

- **cov_param_hat**: Numeric matrix of estimated spatial covariance
  parameters (Nugget, Sill, Range, Smoothness) across samples (*spNNGP
  models only*).

- **phi_hat**: Numeric vector of estimated spatial random effects
  \\\hat{\boldsymbol{\phi}} = \hat{\boldsymbol{\eta}} -
  \boldsymbol{X}\hat{\boldsymbol{\beta}}\\ across all observations.

- **eta_hat**: Numeric vector of latent linear predictors
  \\\hat{\boldsymbol{\eta}}\\ (E-step conditional expectation means).

- **theta_hat**: Numeric vector of estimated Poisson relative rate
  parameters \\\exp(\boldsymbol{X}\hat{\boldsymbol{\beta}} +
  \hat{\boldsymbol{\phi}})\\.

- **predictions**: Numeric vector of predicted counts \\\hat{z} =
  \hat{\theta} \times \text{offset}\\ across all observations.

- **residuals**: Numeric vector of raw residuals (\\z - \hat{z}\\)
  across all observations.

- **beta_tracker**: Numeric matrix recording the trajectory of
  \\\boldsymbol{\beta}\\ estimates across ECM iterations (parameters
  \\\times\\ iterations).

- **gamma_tracker** / **cov_param_tracker**: History of spatial
  parameter estimates across iterations.

- **tau2_tracker**: Numeric matrix recording the trajectory of
  \\\tau^2\\ estimates across iterations (*lattice models only*).

- **R2_tracker**: Numeric matrix of sample-wise \\R^2\\ correlations
  across iterations.

- **MSE_tracker**: Numeric matrix of sample-wise mean squared errors
  across iterations.

- **data_log_like_tracker**: Numeric matrix of observed Poisson data
  log-likelihoods across iterations.

- **expected_log_like_tracker**: Numeric matrix of expected
  complete-data log-likelihoods (\\Q\\-function) across iterations.

- **resid_moran**: Array/matrix recording residual spatial
  autocorrelation (Moran's \\I\\, expected value, and \\p\\-value)
  across iterations.

- **resid_moran_nb**: Matrix recording neighborhood-based residual
  Moran's \\I\\ across iterations (*lattice models only*).

- **beta_neghessian**: Negative Hessian matrix evaluated at
  \\\hat{\boldsymbol{\beta}}\\, representing the observed Fisher
  information used to compute standard errors and Wald statistics.

- **tau2_neghessian**: Vector of negative Hessians for \\\tau^2\\ per
  sample (*lattice models only*).

- **gamma_neghessian**: Vector of negative Hessians for \\\gamma\\ per
  sample (*lattice models only*).

- **performanceSummary**: Data frame summarizing sample-wise model
  performance, goodness-of-fit (normalized MSE, Moran's \\I\\), and
  parameter estimates.

- **start_idx_list**: Named integer vector indicating the starting index
  of each sample in the concatenated observation vectors.

- **time**: `difftime` object recording the total elapsed optimization
  time.

- **run_settings**: List of runtime configuration options, convergence
  thresholds, and initialization flags used during fitting.

## See also

[`TESSERA_lattice`](https://floricaconstantine.github.io/TESSERA/reference/TESSERA_lattice.md),
[`TESSERA_spNNGP`](https://floricaconstantine.github.io/TESSERA/reference/TESSERA_spNNGP.md),
[`calc_Wald_statistics`](https://floricaconstantine.github.io/TESSERA/reference/calc_Wald_statistics.md),
[`summarize_TESSERA`](https://floricaconstantine.github.io/TESSERA/reference/summarize_TESSERA.md)
