# This script aggregates the results of the TESSERA algorithm on the top
# 3000 genes from the data used in the Vignette.
# The script vignette_load_process.R should be run first, before running
# vignette_run_all.R, after which this script can be run.

## Libraries

# General (Numerical) Utilties
library(dplyr)
library(Matrix)


## Set paths

# Paths to data
out_path <- file.path("..", "extdata")
dir.create(out_path, showWarnings = FALSE)
# Edit this
path <- file.path("/scratch", "users", "spatialseq", "natgen_kidney")
in_path <- file.path(path, "vignette_output_data")


## Load data

# List of data outputs from running TESSERA
# Edit if needed
file_list <- list.files(in_path, full.names = TRUE)

performance_df <- list()
wald_df <- list()
for (fname in file_list) {
  # Load in one file
  out <- readRDS(fname)
  
  # Add to lists
  performance_df[[1 + length(performance_df)]] <- out$TESSERA_out$performanceSummary
  out$wald_df$contrast_name <- rownames(out$wald_df)
  wald_df[[1 + length(wald_df)]] <- out$wald_df
}
performance_df <- dplyr::bind_rows(performance_df)
wald_df <- dplyr::bind_rows(wald_df)


## Save processed outputs

saveRDS(
  list(performance_df = performance_df, wald_df = wald_df),
  file.path(out_path, "vignette_processed_output_data2.rds")
)
