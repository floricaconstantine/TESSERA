# Prepare data for the TESSERA method.

Prepare data for the TESSERA method.

## Usage

``` r
prepData(x, ...)

# S4 method for class 'ANY'
prepData(
  x,
  meta_data,
  sample_col,
  design_formula = NULL,
  design_mat = NULL,
  coord_data = NULL,
  D_THRESH = NULL,
  expected_num_neighbors = 6,
  k_search = 20,
  adj_mat = NULL,
  model_type = "Leroux"
)

# S4 method for class 'SpatialExperiment'
prepData(
  x,
  sample_col,
  design_formula = NULL,
  design_mat = NULL,
  model_type = "Leroux",
  D_THRESH = NULL,
  expected_num_neighbors = 6,
  k_search = 20,
  adj_mat = NULL
)
```

## Arguments

- x:

  Variables x measurements (genes x cells) count matrix (can be sparse),
  OR a `SpatialExperiment` object.

- ...:

  Additional arguments passed to methods.

- meta_data:

  Dataframe with metadata/covariates. Ensure that rownames(meta_data) ==
  colnames(x) (if x is a matrix).

- sample_col:

  String: Column name in meta_data identifying which rows correspond to
  which sample (how to break up the data into multiple samples).

- design_formula:

  Formula object (compatible with stats::model.matrix). Supply one of
  design_formula and design_mat.

- design_mat:

  Design matrix. Supply one of design_formula and design_mat. Zero
  columns will be dropped.

- coord_data:

  Dataframe or matrix with coordinates for observations. Rows are
  samples. Ensure that rownames(coord_data) == colnames(x). Needed to
  compute adjacency matrices.

- D_THRESH:

  Distance threshold for assigning two observations as adjacent. Needed
  to compute adjacency matrices.

- expected_num_neighbors:

  In the absence of a known D_THRESH, how many neighbors (on average)
  each cell should have. E.g., Visium data has 6. Ignored if D_THRESH is
  supplied.

- k_search:

  When forming adjacency matrices, the maximum number of neighbors.

- adj_mat:

  A pre-computed adjacency matrix (sparse).

- model_type:

  One of "Leroux", "CAR", "SAR", "spNNGP", or "ALL". Which model will be
  fit by TESSERA. The value "ALL" will prepare the data for all four
  methods.

## Value

A list comprised of the following:

- **coords_list**: A list with a matrix of coordinates for each sample.

- **covariates_list**: A list with a dataframe of metadata for each
  sample.

- **library_size_list**: A list with a vector of library sizes for each
  sample.

- **X_list**: A list with a matrix of design matrices for each sample.

- **W_list**: A list with measurement adjacency matrix for each sample.

- **D_list**: A list with measurement degree matrix for each sample.

- **eig_CS_list**: A list with CAR/SAR model eigenvalues for each
  sample.

- **eig_L_list**: A list with Leroux model eigenvalues for each sample.

## Note

This function supports multiple dispatch for base matrices and
SpatialExperiment objects.

SpatialExperiment support requires: "SpatialExperiment",
"SingleCellExperiment", "SummarizedExperiment".

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
#> Starting Leroux eigenvalue computation for area 1 at 2026-04-15 02:03:38.447104 
#> Starting Leroux eigenvalue computation for area 2 at 2026-04-15 02:03:38.449015 
#> Starting Leroux eigenvalue computation for area 3 at 2026-04-15 02:03:38.453657 

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
