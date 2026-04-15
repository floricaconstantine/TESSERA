# Generate synthetic multi-sample spatial count data

Replaces the observed counts in a `TESSERAData_obj` with synthetic
counts generated from specified true parameters. This utility is
designed for benchmarking and power analyses, allowing users to simulate
data over existing tissue architectures and experimental designs.

## Usage

``` r
prepSynthData(
  TESSERAData_obj,
  gene_list,
  data_gen_model,
  tau2_true = NULL,
  gamma_true = NULL,
  cov_params = NULL,
  cov_type = NULL,
  nngp_k = NULL,
  beta_true = NULL
)
```

## Arguments

- TESSERAData_obj:

  Object containing the experimental structure, typically created by the
  `prepData` method.

- gene_list:

  Character vector: The names of the genes/measurements to simulate.
  These must correspond to row names in the original count data.

- data_gen_model:

  Character: The spatial model to use for simulations. Options include
  "CAR", "SAR", "Leroux" (Lattice models), or "spNNGP".

- tau2_true:

  Spatial scale parameters \\\tau^2\\. Provide a numeric vector (for a
  single gene) or a matrix of dimensions (genes x samples) (for multiple
  genes). Required only for Lattice models. **Note:** \\\tau^2 \ge 0\\
  is required.

- gamma_true:

  Spatial correlation parameters \\\gamma\\. Provide a numeric vector
  (single gene) or a matrix of dimensions (genes x samples). Required
  only for Lattice models.

  - **CAR/SAR**: \\-1 \< \gamma \< 1\\ required.

  - **Leroux**: \\0 \le \gamma \< 1\\ required.

- cov_params:

  Spatial covariance parameters for spNNGP. Provide a matrix (samples x
  3/4 parameters) for a single gene, or an array (genes x samples x 3/4
  parameters) for multiple genes. Parameters must be ordered as (nugget,
  sill, range, and optional smoothness).

- cov_type:

  Character: The spatial correlation kernel for spNNGP. Options are
  "Exp", "Mat", "Gau", and "Sph".

- nngp_k:

  Integer: The number of nearest neighbors to use for the spNNGP kernel
  approximation.

- beta_true:

  True fixed effect coefficients \\\beta\\. Provide a numeric vector
  (single gene) or a matrix of dimensions (genes x covariates).

## Value

A list containing the following components:

- **new_TESSERAData_obj**: A `TESSERAData_obj` where the `counts_list`
  has been replaced with synthetic values.

- **synthetic_count_summary**: A data frame summarizing the generated
  counts and computing similarity metrics between the synthetic and
  original data.

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
#> Starting Leroux eigenvalue computation for area 1 at 2026-04-15 02:03:40.080879 
#> Starting Leroux eigenvalue computation for area 2 at 2026-04-15 02:03:40.082835 
#> Starting Leroux eigenvalue computation for area 3 at 2026-04-15 02:03:40.087712 

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
