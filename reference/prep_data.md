# Prepare data for the TESSERA method.

Prepare data for the TESSERA method.

## Usage

``` r
prep_data(x, ...)

# S4 method for class 'ANY'
prep_data(
  x,
  meta_data,
  sample_col = "sample_id",
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
prep_data(
  x,
  sample_col = "sample_id",
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
# Locate the raw data in inst/extdata
rds_path <- system.file("extdata", "example_raw_data.rds", package = "TESSERA")
# Load the raw data list
raw_data <- readRDS(rds_path)

# Prepare data
TESSERA_data <- prep_data(
  x = raw_data$count_matrix,
  meta_data = raw_data$meta_data,
  sample_col = "sample",
  design_mat = raw_data$design_mat,
  coord_data = raw_data$coords,
  adj_mat = raw_data$W,
  model_type = "Leroux"
)
#> Model type: Leroux
#> Using supplied adjacency matrix.
#> Estimating distance threshold: 0.0904206905630113
#> Subsetting provided adjacency matrix.
#> Starting Leroux eigenvalue computation for area 1 at 2026-08-24 21:11:53.881433 
#> Starting Leroux eigenvalue computation for area 2 at 2026-08-24 21:11:53.884442 
#> Starting Leroux eigenvalue computation for area 3 at 2026-08-24 21:11:53.890227 
```
