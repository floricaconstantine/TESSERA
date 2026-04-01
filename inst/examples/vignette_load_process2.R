# This script processes the data used in the Vignette. 

## Libraries

# TESSERA
library(TESSERA)

# Data Processing and Objects
# BiocManager::install("SummarizedExperiment")
# BiocManager::install("SingleCellExperiment")
# BiocManager::install("SpatialExperiment")
# BiocManager::install("ExperimentHub")
# BiocManager::install("spatialLIBD")
library(SummarizedExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)
library(ExperimentHub)
library(spatialLIBD)

# Data normalization
# BiocManager::install("scran")
library(scran)

# General (Numerical) Utilties
library(dplyr)
library(reshape2)
library(Matrix)
library(MatrixGenerics)


## Load and process data

# Connect to ExperimentHub
ehub <- ExperimentHub::ExperimentHub()
# Download the full real data:
spe_all <- spatialLIBD::fetch_data(type = "spe", eh = ehub)

# Assign row names to be gene names
rownames(spe_all) <- make.unique(SummarizedExperiment::rowData(spe_all)$gene_name, sep="_")

# Add the spatial coordinates into the metadata (convenience)
SummarizedExperiment::colData(spe_all) <- cbind(SummarizedExperiment::colData(spe_all), 
                                                SpatialExperiment::spatialCoords(spe_all))
# Make sure rownames match and are unique
rownames(SummarizedExperiment::colData(spe_all)) <- SummarizedExperiment::colData(spe_all)$key
rownames(SpatialExperiment::spatialCoords(spe_all)) <- SummarizedExperiment::colData(spe_all)$key
colnames(SingleCellExperiment::counts(spe_all)) <- SummarizedExperiment::colData(spe_all)$key

# Specify column names of spatial clusters in colData(spe) 
cluster_col <- "layer_guess_reordered"
# Remove spots missing annotations
spe_all <- spe_all[, !is.na(spe_all[[cluster_col]])]
# Subset to one subject
spe <- spe_all[, colData(spe_all)$subject == "Br8100"]

# Identify mitochondrial genes (using the 'gene_name' column)
is_mito <- grepl("^MT-", rowData(spe)$gene_name, ignore.case = TRUE)
# Subset the object to keep only non-mitochondrial genes
spe <- spe[!is_mito, ]

# Fit the mean-variance trend, blocking by sample_id
dec <- scran::modelGeneVar(spe, block = spe$sample_id)
# Select the top based on biological variance
top_hvgs <- scran::getTopHVGs(dec, n = dim(spe)[1])
# Order spe by HVGs (biological variance)
hvg_order <- order(dec$bio, decreasing = TRUE)
spe <- spe[hvg_order, ]


## Prepare data

TESSERA_data <- prepDataSpatialExperiment(spDataObject = spe, 
                                          sample_col = "sample_id",
                                          design_formula = as.formula( ~ layer_guess_reordered + sample_id),
                                          model_type = "Leroux",
                                          expected_num_neighbors = 6
)


## Save TESSERA_data

out_path <- file.path("..", "extdata")
dir.create(out_path, showWarnings = FALSE)
saveRDS(TESSERA_data,
        file.path(out_path, "vignette_processed_input_data2.rds"))
