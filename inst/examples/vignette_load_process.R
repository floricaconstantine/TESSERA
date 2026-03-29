# This script processes the data used in the Vignette. 

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


## Load and process data

# Load Data
spe <- WeberDivechaLCdata::WeberDivechaLCdata_Visium()
# Keep only spots over tissue
spe <- spe[, SummarizedExperiment::colData(spe)$in_tissue == 1]

# Add the spatial coordinates into the metadata (convenience)
SummarizedExperiment::colData(spe) <- cbind(SummarizedExperiment::colData(spe),
                                            SpatialExperiment::spatialCoords(spe))
# Make sure rownames match
rownames(SpatialExperiment::spatialCoords(spe)) <- rownames(SummarizedExperiment::colData(spe))

# Assign the readable gene names to the rows
rownames(spe) <- SummarizedExperiment::rowData(spe)$gene_name

# Calculate variance across all spots for each gene
rv <- MatrixGenerics::rowVars(SingleCellExperiment::counts(spe))
names(rv) <- rownames(spe)
# Order the entire object by raw count variance (highest to lowest)
spe <- spe[order(rv, decreasing = TRUE), ]


## Prepare Data

# Prepare data for the TESSERA package
# Simple design: Look at annotated region in v. out
TESSERA_data <- TESSERA::prepDataSpatialExperiment(
  spDataObject = spe,
  sample_col = "sample_part_id",
  design_formula = as.formula( ~ annot_region + sample_part_id),
  model_type = "Leroux",
  expected_num_neighbors = 6
)


## Save TESSERA_data

out_path <- file.path("..", "extdata")
dir.create(out_path, showWarnings = FALSE)
saveRDS(TESSERA_data,
        file.path(out_path, "vignette_processed_input_data.rds"))
