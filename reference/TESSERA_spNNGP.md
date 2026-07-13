# Fit Multi-Sample Poisson Spatial GLMM via spNNGP

Fits a multi-sample Poisson spatial generalized linear mixed model
(GLMM) using a shared set of fixed effects across all samples while
permitting sample-specific spatial random effects modeled via a Sparse
Nearest Neighbor Gaussian Process (spNNGP). Supports Exponential,
Matern, Gaussian, and Spherical covariance kernels.

## Usage

``` r
TESSERA_spNNGP(
  TESSERAData_obj,
  gene_name,
  cov_type = "Exp",
  nngp_k = 20,
  em_iters = 200,
  opt_iters = 5,
  em_min_iters = 15,
  em_tol = 0.001,
  em_stopping = NULL,
  beta_init = "glm",
  cov_init = "BRISC",
  cov_fit_method = "BRISC",
  verbose = FALSE,
  dense_matrices = FALSE
)
```

## Arguments

- TESSERAData_obj:

  An object containing prepared data, typically created by
  [`prep_data`](https://floricaconstantine.github.io/TESSERA/reference/prep_data.md).

- gene_name:

  Character: The name of the gene/measurement (row) to fit.

- cov_type:

  Character: The covariance kernel for the Gaussian process. Options are
  "Exp", "Mat", "Gau", and "Sph".

- nngp_k:

  Integer: The number of nearest neighbors to use for the spNNGP kernel
  approximation.

- em_iters:

  Integer: Maximum number of ECM iterations.

- opt_iters:

  Integer: Number of inner CM steps (Conditional Maximization) per
  iteration. The algorithm maximizes the expected likelihood for the
  covariance parameters, then \\\beta\\, holding other parameters
  constant.

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

- cov_init:

  Method for initializing covariance parameters. Options: "BRISC" (calls
  BRISC on residuals), "variogram" (initializes via empirical
  variogram), or a numeric vector/matrix.

- cov_fit_method:

  Character: The method used for the M-step update of covariance
  parameters. Options are "BRISC" or "variogram".

- verbose:

  Logical: Whether to print iteration-wise parameter updates.

- dense_matrices:

  Logical: If `TRUE`, treats the precision matrix \\Q\\ as dense during
  specific E-step calculations. This dramatically increases memory usage
  but can lead to a 60 to 80 percent decrease in computation time.

## Value

A list containing the following components:

- **beta_hat**: Estimated fixed effect coefficients.

- **cov_param_hat**: Estimated spatial covariance parameters (nugget,
  sill, range, and optionally smoothness).

- **phi_hat**: Estimated spatial random effects.

- **theta_hat**: Estimated Poisson rate parameters \\exp(X\beta +
  \phi)\\.

- **beta_tracker**: Matrix of \\\beta\\ estimates across iterations.

- **cov_param_tracker**: Array of covariance parameter history (areas x
  iterations x parameters).

- **data_log_like_tracker**: Total log-likelihood across iterations.

- **MSE_tracker**: Mean Squared Error history across iterations.

- **beta_neghessian**: Negative Hessian of the final \\\beta\\ estimate.

- **performanceSummary**: A summary data frame for each sample; see
  [`summarize_TESSERA`](https://floricaconstantine.github.io/TESSERA/reference/summarize_TESSERA.md).

- **time**: Total execution time.

- **run_settings**: List of parameters used for the run.

## Note

The spNNGP approach is designed for scalability in datasets where
traditional lattice adjacency is difficult to define.

## References

Meng, Xiao-Li, and Donald B. Rubin. "Maximum likelihood estimation via
the ECM algorithm: A general framework." Biometrika 80.2 (1993):
267-278.

Saha, Arkajyoti, and Abhirup Datta. "BRISC: bootstrap for rapid
inference on spatial covariances." Stat 7.1 (2018): e184.

## Author

Florica J Constantine, florica AT berkeley.edu

## Examples

``` r
# Locate the prepped TESSERAData object in inst/extdata
rds_path <- system.file("extdata", "example_prepData.rds", package = "TESSERA")
# Load the TESSERAData object
TESSERA_data <- readRDS(rds_path)

# Fit the Poisson generalized spatial linear model using spNNGP
TESSERA_out_spNNGP <- suppressMessages(suppressWarnings(
  TESSERA_spNNGP(
    TESSERAData_obj = TESSERA_data,
    gene_name = "example",
    cov_type = "Exp",
    em_iters = 2,
    opt_iters = 1,
    verbose = FALSE
  )
))
```
