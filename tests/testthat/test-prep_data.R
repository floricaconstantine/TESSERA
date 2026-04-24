library(testthat)
library(Matrix)
library(ggplot2)

# --- Shared Setup for Data Prep Tests ---
set.seed(2026)
n_cells <- 10
genes <- c("GeneA", "GeneB")
counts <- matrix(rpois(n_cells * 2, 5), nrow = 2)
rownames(counts) <- genes
colnames(counts) <- paste0("Cell", 1:10)

meta <- data.frame(
  sample = rep(c("S1", "S2"), each = 5),
  batch = rep(1:2, 5),
  row.names = colnames(counts)
)

coords <- matrix(runif(n_cells * 2), ncol = 2)
rownames(coords) <- colnames(counts)

# --- prep_data Tests ---

test_that("prep_data validates design and dimension arguments", {
  # Must supply exactly one of design_formula or design_mat
  expect_error(prep_data(counts, meta, "sample"),
               "Must supply one of design formula")
  expect_error(
    prep_data(
      counts,
      meta,
      "sample",
      design_formula = ~ batch,
      design_mat = matrix(1, 10, 1)
    ),
    "Cannot supply both design formula and design matrix"
  )
  
  # Dimension mismatch (counts vs meta)
  expect_error(prep_data(counts[, 1:5], meta, "sample", design_formula = ~ batch))
})

test_that("prep_data creates a valid TESSERAData object with multiple samples",
          {
            # suppressMessages silences adjacency and eigenvalue logs
            res <- suppressMessages(
              prep_data(
                x = counts,
                meta_data = meta,
                sample_col = "sample",
                design_formula = ~ batch,
                coord_data = coords,
                D_THRESH = 0.5,
                model_type = "Leroux"
              )
            )
            
            expect_s3_class(res, "TESSERAData")
            expect_length(res$counts_list, 2)
            expect_named(res$counts_list, c("S1", "S2"))
            
            # Verify Leroux eigenvalues were computed
            expect_false(is.null(res$eig_L_list))
            expect_length(res$eig_L_list$S1, 5)
          })

test_that("prep_data handles model_type dispatch and match.arg logic", {
  # Test "ALL" model type
  res_all <- suppressMessages(
    prep_data(
      counts,
      meta,
      "sample",
      design_formula = ~ 1,
      coord_data = coords,
      D_THRESH = 0.5,
      model_type = "ALL"
    )
  )
  expect_false(is.null(res_all$eig_L_list))
  expect_false(is.null(res_all$eig_CS_list))
  
  # Test "spNNGP" which should set compute_eigs to "NONE"
  expect_warning(suppressMessages(
    prep_data(
      counts,
      meta,
      "sample",
      design_formula = ~ 1,
      coord_data = coords,
      D_THRESH = 0.5,
      model_type = "spNNGP"
    )
  ),
  "Not computing eigenvalues: option selected")
  
  # Test invalid model type triggers match.arg error
  expect_error(
    prep_data(
      counts,
      meta,
      "sample",
      design_formula = ~ 1,
      model_type = "InvalidType"
    ),
    "'arg' should be one of"
  )
})

test_that("prep_data automatically estimates D_THRESH if missing", {
  # We use expected_num_neighbors = 2 because each sample only has 5 cells
  expect_message(suppressWarnings(
    prep_data(
      counts,
      meta,
      "sample",
      design_formula = ~ 1,
      coord_data = coords,
      expected_num_neighbors = 2
    )
  ),
  "Estimating distance threshold")
})

test_that("prep_data handles library size normalization logic", {
  # Pass coord_data and D_THRESH to avoid the automatic thresholding crash
  res_multi <- suppressWarnings(suppressMessages(
    prep_data(
      counts,
      meta,
      "sample",
      design_formula = ~ 1,
      coord_data = coords,
      D_THRESH = 0.5
    )
  ))
  expect_equal(as.numeric(res_multi$library_size_list$S1[1]), sum(counts[, 1]))
  
  # Single-gene: library size is set to 1
  res_single <- suppressWarnings(suppressMessages(
    prep_data(
      counts[1, , drop = FALSE],
      meta,
      "sample",
      design_formula = ~ 1,
      coord_data = coords,
      D_THRESH = 0.5
    )
  ))
  expect_equal(as.numeric(res_single$library_size_list$S1[1]), 1)
})

test_that("prep_data triggers warnings for poor design matrices", {
  bad_mat_zero <- matrix(c(rep(1, 5), rep(0, 5)), ncol = 2)
  rownames(bad_mat_zero) <- colnames(counts)[1:5]
  
  expect_warning(
    expect_warning(suppressMessages(
      prep_data(
        counts[, 1:5],
        meta[1:5, ],
        "sample",
        design_mat = bad_mat_zero,
        coord_data = coords[1:5, ],
        D_THRESH = 0.5
      )
    ), "all-zero columns", fixed = TRUE),
    "Dropping columns",
    fixed = TRUE
  )
  
  bad_mat_rank <- matrix(rep(1, 10), ncol = 2)
  rownames(bad_mat_rank) <- colnames(counts)[1:5]
  
  expect_warning(suppressMessages(
    prep_data(
      counts[, 1:5],
      meta[1:5, ],
      "sample",
      design_mat = bad_mat_rank,
      coord_data = coords[1:5, ],
      D_THRESH = 0.5
    )
  ), "not full rank", fixed = TRUE)
})

# --- plot_neighbor_distances Tests ---

test_that("plot_neighbor_distances returns the correct visualization list",
          {
            # Use k_search = 3 because each sample only has 5 cells
            res <- plot_neighbor_distances(meta, "sample", coords, k_search = 3)
            
            expect_named(res,
                         c("nb_dist", "mean_nb_dist", "std_nb_dist", "distances_plot"))
            expect_s3_class(res$distances_plot, "ggplot")
            expect_equal(nrow(res$mean_nb_dist), 2)
          })

# --- SpatialExperiment Method Tests ---

test_that("prep_data handles SpatialExperiment objects via S4 dispatch", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("SingleCellExperiment")
  
  # Create mock SpatialExperiment
  se <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    colData = meta,
    spatialCoords = coords
  )
  
  # Dispatch to SpatialExperiment method
  # Use expected_num_neighbors = 2 to avoid nndist extent errors
  res <- suppressMessages(
    prep_data(
      x = se,
      sample_col = "sample",
      design_formula = ~ batch,
      model_type = "Leroux",
      D_THRESH = 0.5,
      expected_num_neighbors = 2
    )
  )
  
  expect_s3_class(res, "TESSERAData")
  expect_length(res$counts_list, 2)
  expect_false(is.null(res$eig_L_list))
})
