# This script runs the TESSERA algorithm on the top 3000 genes from the data 
# used in the Vignette. 
# It is written to be run on an HPC computing cluster.
# The script vignette_load_process.R should be run first. 

## Libraries

# TESSERA
library(TESSERA)

# General (Numerical) Utilties
library(dplyr)
library(Matrix)


## Load in data

# Paths to data
in_path <- file.path("..", "extdata")
# Edit this
path <- file.path("/scratch", "users", "spatialseq", "natgen_kidney")
out_path <- file.path(path, "vignette_output_data")

# Load prepData object created in other script (vignette_load_process.R)
TESSERA_data <- readRDS(file.path(in_path, "vignette_processed_input_data2.rds"))


## Parameters

em_iters <- 200
opt_iters <- 5
em_min_iters <- 15
em_tol <- 1e-3
em_stopping <- "rel_loglike"
em_verbose <- TRUE 
cov_type <- "Mat"

# Key markers for norepinephrine neurons
gene_names <- rownames(TESSERA_data$counts_list[[1]])

# Contrast matrix
layer_pairs <- combn(7, 2)
contrast_matrix <- matrix(
  data = 0,
  ncol = ncol(TESSERA_data$X_list[[1]]),
  nrow = ncol(layer_pairs)
)

rownames_contrast <- rep("", ncol(layer_pairs))
for (idx in 1:ncol(layer_pairs)) {
  # This model has an intercept, so all effects are relative to layer 1
  if (1 == min(layer_pairs[, idx])) {
    contrast_matrix[idx, max(layer_pairs[, idx])] <- 1
  } else {
    contrast_matrix[idx, max(layer_pairs[, idx])] <- 1
    contrast_matrix[idx, min(layer_pairs[, idx])] <- -1
  }
  rownames_contrast[idx] <- paste0("L", max(layer_pairs[, idx]), "-", "L", min(layer_pairs[, idx]))
}
rownames(contrast_matrix) <- rownames_contrast
colnames(contrast_matrix) <- colnames(TESSERA_data$X_list[[1]])

## Run Settings

# Which models to fit
fit_model_list <- c("Leroux")
# Which genes to fit on
gene_index_list <- seq_len(5000)
# Data frame with run settings
sim_df <- expand.grid(gene_index = gene_index_list, fit_model = fit_model_list)

# Read in index of simulation job from command line input
run_idx <- NA
offset <- 0
args <- commandArgs(trailingOnly = TRUE)
if (0 == length(args)) {
  stop("No arguments supplied: ERROR")
} else {
  for (idx in 1:length(args)) {
    eval(parse(text = args[idx]))
  }
  
  print(run_idx)
  print(offset)
  run_idx <- run_idx + offset
  print(run_idx)
  
  # run_idx must be defined in arguments
  if ((1 > run_idx) | (nrow(sim_df) < run_idx)) {
    stop("run_idx must be between 1 and # of possible jobs")
  }
}

# Select model and gene
fit_model <- sim_df[run_idx, ]$fit_model
gene_idx <- sim_df[run_idx, ]$gene_index
cat(fit_model, gene_idx, gene_names[gene_idx], "\n")

## Run TESSERA

if ("Leroux" == fit_model) {
  TESSERA_out <- TESSERA::TESSERA_lattice(
    TESSERAData_obj = TESSERA_data,
    gene_name = gene_names[gene_idx],
    model_type = "Leroux",
    em_iters = em_iters,
    # How many ECM iterations to run
    opt_iters = opt_iters,
    # How many CM steps to run per E-Step
    em_min_iters = em_min_iters,
    # Min ECM iterations before early stopping
    em_tol = em_tol,
    # Tolerance for early stopping
    em_stopping = em_stopping,
    # How to determine early stopping: Relative change in log likelihood
    verbose = em_verbose
  )
} else if ("spNNGP" == fit_model) {
  TESSERA_out <- TESSERA::TESSERA_spNNGP(
    TESSERAData_obj = TESSERA_data,
    gene_name = gene_names[gene_idx],
    cov_type = cov_type,
    em_iters = em_iters,
    # How many ECM iterations to run
    opt_iters = opt_iters,
    # How many CM steps to run per E-Step
    em_min_iters = em_min_iters,
    # Min ECM iterations before early stopping
    em_tol = em_tol,
    # Tolerance for early stopping
    em_stopping = em_stopping,
    # How to determine early stopping: Relative change in log likelihood
    verbose = em_verbose
  )
}


## Inference

wald_df <- TESSERA::waldTestStastics(TESSERA_out, contrast_matrix)
# Drop column for easy viewing
wald_df <- dplyr::select(wald_df, -c("contrast_string"))

## Save

saveRDS(list(TESSERA_out = TESSERA_out, wald_df = wald_df),
        file.path(
          out_path,
          paste0(
            "vignette_raw_output_",
            "gene_",
            gene_names[gene_idx],
            "_fitmodel_",
            fit_model,
            ".rds"
          )
        ))
