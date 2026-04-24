# This script processes the data used in the Vignette.

## Libraries

# TESSERA
library(TESSERA)

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


## Load data

# The data is already in a SpatialExperiment object.
# The rows (counts) are sorted by variance (decreasing).
# The rownames of the metadata and of the coordinates match the
# colnames of the counts.
spe <- readRDS(file.path(".", "kidney_data_spe.rds"))


## Create design matrix

# We use a design with no intercept and effects for all pairs of interactions of
# cell type and condition.
# We also include effects for sample; this causes some complexity, as shown below.

# Which covariates to use
celltype_col <- "celltype"
sample_col <- "orig.ident"
group_col <- "Group"
coord_cols <- c("pxl_row_in_fullres", "pxl_col_in_fullres")
covariates_keep <- c(group_col, celltype_col) # Use when creating design formula

# Reference celltype, sample, and condition
ref_celltype <- "C_TAL"
ref_sample <- "HK2874_ST"
ref_group <- "Control"

# Re-level factors
SummarizedExperiment::colData(spe)[, sample_col] <- relevel(factor(SummarizedExperiment::colData(spe)[, sample_col]), ref = ref_sample)
SummarizedExperiment::colData(spe)[, celltype_col] <- relevel(factor(SummarizedExperiment::colData(spe)[, celltype_col]), ref = ref_celltype)
SummarizedExperiment::colData(spe)[, group_col] <- relevel(factor(SummarizedExperiment::colData(spe)[, group_col]), ref = ref_group)

## If we want a sample/patient effect; we need to nest this inside disease
# Create a new 'nested' variable as in the DESeq2 vignette
control_ids <- unique(SummarizedExperiment::colData(spe)[, sample_col][SummarizedExperiment::colData(spe)[, group_col]  == "Control"])
DKD_ids <- unique(SummarizedExperiment::colData(spe)[, sample_col][SummarizedExperiment::colData(spe)[, group_col]  == "DKD"])
HKD_ids <- unique(SummarizedExperiment::colData(spe)[, sample_col][SummarizedExperiment::colData(spe)[, group_col]  == "HKD"])
new_id <- rep(0, nrow(SummarizedExperiment::colData(spe)))
for (idx in 1:length(control_ids)) {
  new_id[SummarizedExperiment::colData(spe)[, sample_col] == control_ids[idx]] <- idx
}
for (idx in 1:length(DKD_ids)) {
  new_id[SummarizedExperiment::colData(spe)[, sample_col] == DKD_ids[idx]] <- idx
}
for (idx in 1:length(HKD_ids)) {
  new_id[SummarizedExperiment::colData(spe)[, sample_col] == HKD_ids[idx]] <- idx
}
SummarizedExperiment::colData(spe)$nested_id <- factor(new_id)
new_nested_sample_col <- "nested_id"

# Create dataframe
data_mat <- SummarizedExperiment::colData(spe)
data_mat$z <- SingleCellExperiment::counts(spe)[1, ]

# Design formula
design_formula <- as.formula(paste0("z ~ ", "0 + ", paste(
  covariates_keep[1],
  c(covariates_keep[-c(1)], new_nested_sample_col),
  sep = ":",
  collapse = " + "
)))
print(design_formula)
# Create design matrix
X_mat <- model.matrix(design_formula, data_mat)

# Rename sample/patient effects for clarity later
for (idx in 2:length(DKD_ids)) {
  colnames(X_mat) <- gsub(
    paste0("GroupDKD:nested_id", idx),
    paste0("GroupDKD:", DKD_ids[idx]),
    colnames(X_mat)
  )
}
for (idx in 2:length(HKD_ids)) {
  colnames(X_mat) <- gsub(
    paste0("GroupHKD:nested_id", idx),
    paste0("GroupHKD:", HKD_ids[idx]),
    colnames(X_mat)
  )
}
for (idx in 2:length(control_ids)) {
  colnames(X_mat) <- gsub(
    paste0("GroupControl:nested_id", idx),
    paste0("GroupControl:", control_ids[idx]),
    colnames(X_mat)
  )
}

# Toss out zero columns
X_mat <- X_mat[, colSums(X_mat) != 0]
# Make sure design is full rank
stopifnot(Matrix::rankMatrix(X_mat)[1] == ncol(X_mat))

## Prepare Data

# Prepare data for the TESSERA package
# Simple design: Look at annotated region in v. out
TESSERA_data <- TESSERA::prep_data(
  x = spe,
  sample_col = sample_col,
  design_mat = X_mat,
  model_type = "Leroux",
  expected_num_neighbors = 6
)


## Save prepared data object

# Use "xz" compression for smallest file size
saveRDS(TESSERA_data,
        file.path(".", "kidney_data_processed_TESSERA.rds"),
        compress = "xz")
