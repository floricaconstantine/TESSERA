# Fit Multi-Sample Poisson Spatial GLMM via ECM Algorithm

Fits a multi-sample Poisson spatial generalized linear mixed model
(GLMM) using a shared set of fixed effects across all samples while
permitting sample-specific spatial random effects (CAR, SAR, or Leroux).
Parameter estimation is performed via an Expectation-Conditional
Maximization (ECM) algorithm.

## Usage

``` r
TESSERA_lattice(
  TESSERAData_obj,
  gene_name,
  model_type = "Leroux",
  em_iters = 200,
  opt_iters = 5,
  em_min_iters = 15,
  em_tol = 0.001,
  em_stopping = NULL,
  beta_init = "glm",
  gamma_init = "moran",
  tau2_init = "var",
  verbose = TRUE,
  dense_matrices = FALSE
)
```

## Arguments

- TESSERAData_obj:

  An object containing prepared data, typically created by
  [`prep_data`](https://floricaconstantine.github.io/TESSERA/reference/prep_data.md).

- gene_name:

  Character: The name of the gene/measurement (row) to fit.

- model_type:

  Character: The spatial model for random effects. Options are "CAR",
  "SAR", or "Leroux".

- em_iters:

  Integer: Maximum number of ECM iterations.

- opt_iters:

  Integer: Number of inner CM steps (Conditional Maximization) per
  iteration. The algorithm maximizes the expected likelihood for
  \\\tau^2\\, then \\\gamma\\ for each area, and finally \\\beta\\ while
  holding other parameters constant.

- em_min_iters:

  Integer: Minimum number of ECM iterations to perform before allowing
  early stopping.

- em_tol:

  Numeric: Convergence tolerance for early stopping.

- em_stopping:

  Character: Metric used for early stopping:

  - `NULL`: No early stopping.

  - "abs_loglike": Absolute change in total log-likelihood.

  - "rel_loglike": Relative change in total log-likelihood.

  - "abs_beta_norm": Absolute \\L_2\\ norm of the change in \\\beta\\.

  - "rel_beta_norm": Relative \\L_2\\ norm of the change in \\\beta\\.

- beta_init:

  Initial value for \\\beta\\. Options: "glm" (fit a Poisson GLM),
  "random" (standard normal), or a numeric vector.

- gamma_init:

  Initial value for \\\gamma\\. Options: "moran" (absolute Moran's I of
  residuals), "random" (standard uniform), or numeric.

- tau2_init:

  Initial value for \\\tau^2\\. Options: "lognormal" (approximate
  Poisson-lognormal variance), "var" (variance of residuals), "random",
  or numeric.

- verbose:

  Logical: Whether to print iteration-wise parameter updates.

- dense_matrices:

  Logical: If `TRUE`, treats the precision matrix \\Q\\ as dense during
  specific E-step calculations. This increases memory usage but can
  improve computation speed by 10 to 20 percent.

## Value

A list containing the following components:

- **beta_hat**: Estimated fixed effect coefficients.

- **gamma_hat**: Estimated spatial correlation parameters (per area).

- **tau2_hat**: Estimated spatial scale parameters (per area).

- **phi_hat**: Estimated spatial random effects.

- **theta_hat**: Estimated Poisson rate parameters \\exp(X\beta +
  \phi)\\.

- **beta_tracker**: Matrix of \\\beta\\ estimates across iterations.

- **gamma_tracker**: Matrix of \\\gamma\\ estimates across iterations.

- **tau2_tracker**: Matrix of \\\tau^2\\ estimates across iterations.

- **performanceSummary**: A summary data frame for each sample.

- **time**: Total execution time.

## Note

For CAR and SAR models, ensure the adjacency structure \\W\\ contains no
isolated points (every observation must have \\\ge 1\\ neighbor) to
ensure invertibility.

## References

Meng, Xiao-Li, and Donald B. Rubin. "Maximum likelihood estimation via
the ECM algorithm: A general framework." Biometrika 80.2 (1993):
267-278.

## Author

Florica J Constantine, florica AT berkeley.edu

## Examples

``` r
# Locate the prepped TESSERAData object in inst/extdata
rds_path <- system.file("extdata", "example_prepData.rds", package = "TESSERA")
# Load the TESSERAData object
TESSERA_data <- readRDS(rds_path)

# Fit the Poisson generalized spatial linear model
TESSERA_out_Leroux <- suppressMessages(suppressWarnings(
  TESSERA_lattice(
    TESSERAData_obj = TESSERA_data,
    gene_name = "example",
    model_type = "Leroux",
    em_iters = 2,
    opt_iters = 1,
    verbose = FALSE
  )
))
```
