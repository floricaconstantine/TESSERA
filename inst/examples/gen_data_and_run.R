## Script to generate data for runnable examples. 

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

saveRDS(
  list(
    count_matrix = count_matrix,
    meta_data = meta_data,
    design_mat = design_mat,
    coords = coords,
    W = W
  ),
  file.path("..", "extdata", "example_raw_data.rds")
)

# Prepare data
TESSERA_data <- prep_data(
  x = count_matrix,
  meta_data = meta_data,
  sample_col = "sample",
  design_mat = design_mat,
  coord_data = coords,
  adj_mat = W,
  model = "Leroux"
)

saveRDS(TESSERA_data,
        file.path("..", "extdata", "example_prepData.rds"))


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

saveRDS(TESSERA_out_Leroux,
        file.path("..", "extdata", "example_TESSERA_out_Leroux.rds"))


## Inference

# To get Wald statistics
calc_Wald_statistics(TESSERAOutput_obj = TESSERA_out_Leroux,
                     contrast_mat = matrix(c(1, 0, 0), nrow = 1))


## Resampling / parametric data generation

# Resample the data using the same coordinates and covariates
# but with the fitted parameters
TESSERA_resampled_data <- prep_synth_data(
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
