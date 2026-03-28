## Libraries

# TESSERA
library(TESSERA)

# Data
# BiocManager::install("WeberDivechaLCdata")
library(WeberDivechaLCdata)

# Data Processing and Objects
# BiocManager::install("SummarizedExperiment")
# BiocManager::install("SingleCellExperiment")
# BiocManager::install("SpatialExperiment")
library(SummarizedExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)

# General (Numerical) Utilties
library(dplyr)
library(reshape2)
library(Matrix)
library(MatrixGenerics)
library(MASS)


## Load in data

# Paths to data
path <- file.path("/scratch", "users", "spatialseq", "natgen_kidney")
in_path <- file.path(path, "vignette_input_data")
out_path <- file.path(path, "vignette_output_data")

# Load prepData object created in other script (vignette_load_process.R)
TESSERA_data <- readRDS(file.path(in_path, "vignette_processed_input_data.rds"))


## Parameters

em_iters <- 100
opt_iters <- 5
em_min_iters <- 15
em_tol <- 1e-3
em_stopping <- "rel_loglike"
em_verbose <- FALSE
cov_type <- "Mat"

# Key markers for norepinephrine neurons
gene_names <- rownames(TESSERA_data$counts_list[[1]])
lc_markers <- c("TH", "DBH", "SLC6A2", "SLC18A2")

# Contrast matrix
contrast_matrix <- matrix(
  data = 0,
  ncol = ncol(TESSERA_data$X_list[[1]]),
  nrow = 1
)
colnames(contrast_matrix) <- colnames(TESSERA_data$X_list[[1]])
rownames(contrast_matrix) <- colnames(TESSERA_data$X_list[[1]])[2]
contrast_matrix[1, 2] <- 1


## Run Settings

# Check which marker genes are in the data set
present_markers <- intersect(lc_markers, gene_names)
marker_gene_idx <- which(gene_names %in% lc_markers)

# Which models to fit
fit_model_list <- c("Leroux", "spNNGP")
# Which genes to fit on
gene_index_list <- unique(sort(c(seq_len(3000), marker_gene_idx)))
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
