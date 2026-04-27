# Generate synthetic multi-sample spatial count data

Replaces the observed counts in a `TESSERAData_obj` with synthetic
counts generated from specified true parameters. This utility is
designed for benchmarking and power analyses, allowing users to simulate
data over existing tissue architectures and experimental designs.

## Usage

``` r
prep_synth_data(
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
# Locate the prepped TESSERAData object in inst/extdata
rds_path <- system.file("extdata", "example_prepData.rds", package = "TESSERA")
# Load the TESSERAData object
TESSERA_data <- readRDS(rds_path)

#' # Locate the saved model results in inst/extdata
rds_path <- system.file("extdata", "example_TESSERA_out_Leroux.rds", package = "TESSERA")
# Load the results object
TESSERA_out_Leroux <- readRDS(rds_path)

# Run the synthetic data generation
TESSERA_resampled_data <- prep_synth_data(
  TESSERAData_obj = TESSERA_data,
  gene_list = "example",
  data_gen_model = "Leroux",
  tau2_true = TESSERA_out_Leroux$tau2_hat,
  gamma_true = TESSERA_out_Leroux$gamma_hat,
  beta_true = TESSERA_out_Leroux$beta_hat
)$new_TESSERAData_obj
```
