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
  dense_matrices = TRUE
)
```

## Arguments

- TESSERAData_obj:

  An object containing prepared data, typically created by
  [`prepData`](prepData.md).

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
  [`summarizeTESSERAPerformance`](summarizeTESSERAPerformance.md).

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
# Load package
library(TESSERA)

set.seed(2026)

## Generate synthetic data

data_list <- list()
n_samples <- 3 #
gamma_list <- c(0.1, 0.5, 0.9) # Spatial correlation parameter
tau2_list <- c(2.0, 0.1, 1.0) # Spatial variance parameter
beta_true <- c(1, 0, -1) # True fixed effects
n_list <- c(200, 300, 400) # Number of measurements per sample
nb_dist <- 0.05 # Distance to determine if a point is a neighbor
for (idx in seq_len(n_samples)) {
  data_list[[idx]] <- generate_data_one_area(
    n_points = n_list[idx],
    nb_dist = nb_dist,
    model_type = "Leroux",
    beta_true = beta_true,
    gamma_true = gamma_list[idx],
    tau2_true = tau2_list[idx]
  )
}

## Prepare for TESSERA

# This is how you'd go about preparing data manually
# See the prepData for a much more streamlined workflow with 
# a SpatialExperiment object

# Count matrix: one gene per row
count_matrix <- matrix(Reduce(c, sapply(data_list, function (x) {
  x$z
})), nrow = 1)
rownames(count_matrix) <- c("example")
colnames(count_matrix) <- 1:ncol(count_matrix)
# Meta data
meta_data <- data.frame(sample = Reduce(c, sapply(seq_len(n_samples), function (x) {
  rep(x, n_list[x])
})))
rownames(meta_data) <- colnames(count_matrix)
# Coordinates
coords <- Reduce(rbind, sapply(data_list, function (x) {
  x$coords
}))
rownames(coords) <- colnames(count_matrix)
# Design matrix
design_mat <- Reduce(rbind, sapply(data_list, function (x) {
  x$X
}))
rownames(design_mat) <- colnames(count_matrix)
# Adjacency matrix
W <- Matrix::bdiag(sapply(data_list, function (x) {
  x$W
}))
rownames(W) <- colnames(count_matrix)
colnames(W) <- colnames(count_matrix)

# Prepare data
TESSERA_data <- prepData(
  x = count_matrix,
  meta_data = meta_data,
  sample_col = "sample",
  design_mat = design_mat,
  coord_data = coords,
  adj_mat = W,
  model = "Leroux"
)
#> Model type: Leroux
#> Using supplied adjacency matrix.
#> Estimating distance threshold: 0.0904206905630113
#> Subsetting provided adjacency matrix.
#> Starting Leroux eigenvalue computation for area 1 at 2026-04-15 02:03:31.373405 
#> Starting Leroux eigenvalue computation for area 2 at 2026-04-15 02:03:31.375286 
#> Starting Leroux eigenvalue computation for area 3 at 2026-04-15 02:03:31.379946 

# Optional checks for data objects
checkInputsTESSERA(TESSERA_data)
checkInputsTESSERAspNNGP(TESSERA_data)

## Fit models

# Note: In these examples we only run 2 iterations with 1 inner CM iteration.
# In practice, more iterations are needed.
suppressMessages(suppressWarnings({
  # Run TESSERA with a Leroux (lattice) model
  TESSERA_out_Leroux <- TESSERA_lattice(
    TESSERAData_obj = TESSERA_data,
    gene_name = "example",
    model_type = "Leroux",
    em_iters = 2,
    opt_iters = 1,
    verbose = FALSE
  )
  
  # Run TESSERA with a spNNGP (sparse Gaussian kernel) model
  TESSERA_out_spNNGP <- TESSERA_spNNGP(
    TESSERAData_obj = TESSERA_data,
    gene_name = "example",
    cov_type = "Mat",
    em_iters = 2,
    opt_iters = 1,
    verbose = FALSE
  )
}))

## Inference

# To get Wald statistics
waldTestStastics(TESSERAOutput_obj = TESSERA_out_Leroux,
                 contrast_mat = matrix(c(1, 0, 0), nrow = 1))
#>      gene fit_model kernel_type contrast_string contrast_indices contrast_val
#> 1 example    Leroux          NA        1*+0*+0*                1     1.021903
#>   contrast_se wald_stat_t
#> 1  0.05245812    19.48037

## Resampling / parametric data generation

# Resample the data using the same coordinates and covariates
# but with the fitted parameters
TESSERA_resampled_data <- prepSynthData(
  TESSERAData_obj = TESSERA_data,
  gene_list = "example",
  data_gen_model = "Leroux",
  tau2_true = TESSERA_out_Leroux$tau2_hat,
  gamma_true = TESSERA_out_Leroux$gamma_hat,
  beta_true = TESSERA_out_Leroux$beta_hat
)$new_TESSERAData_obj
# Refit on the synthetically generated data
# Note: In these examples we only run 2 iterations with 1 inner CM iteration.
# In practice, more iterations are needed.
suppressMessages(suppressWarnings({
  TESSERA_out_Leroux <- TESSERA_lattice(
    TESSERAData_obj = TESSERA_resampled_data,
    gene_name = "example",
    model_type = "Leroux",
    em_iters = 2,
    opt_iters = 1,
    verbose = FALSE
  )
}))
```
