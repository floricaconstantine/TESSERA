## Purpose of file: Prepare the data to run the TESSERA algorithm
# Create various lists containing, respectively,
# coordinate matrices, covariate matrices, count matrices, library sizes, etc.
# with a list element for each sample
## Florica Constantine (florica AT berkeley DOT edu)

#' Prepare data for the TESSERA method.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param count_matrix Variables x samples (genes x cells) count matrix.
#'  Can be sparse. Rownames are variable names.
#'  Ensure that there is at least one row (a matrix).
#' @param meta_data Dataframe with metadata/covariates.
#'  Ensure that rownames(meta_data) == colnames(count_matrix).
#' @param sample_col String: Column name in meta_data identifying which rows correspond
#'  to which sample (how to break up the data into multiple samples).
#' @param design_formula Formula object (compatible with stats::model.matrix).
#'  Supply one of design_formula and design_mat.
#' @param design_mat Design matrix.
#'  Supply one of design_formula and design_mat.
#'  Zero columns will be dropped.
#' @param coord_data Dataframe or matrix with coordinates for observations.
#'  Rows are samples.
#'  Ensure that rownames(coord_data) == colnames(count_matrix).
#'  Needed to compute adjacency matrices.
#'  Needed for spNNGP model, not necessary for lattice models if adjacency matrices
#'  are supplied elsewhere.
#' @param D_THRESH Distance threshold for assigning two observations as adjacent.
#'  Needed to compute adjacency matrices.
#'  Not needed for non-lattice models (set to NULL).
#' @param k_search When forming adjacency matrices, the maximum number of neighbors.
#' @param adj_mat A pre-computed adjacency matrix (sparse).
#'  Instead of forming it via D_THRESH or coordinates, a user can supply it.
#' @param compute_eigs Whether to compute eigenvalues of various matrices.
#'  CAR/SAR: D^\{-1\} W; Leroux: D - W.
#'  Depends on forming adjacency matrices, so needs coord_data, D_THRESH, k_search.
#'  "None": Compute no eigenvalues.
#'  String containing "C": Compute eigenvalues for the CAR model.
#'  String containing "S": Compute eigenvalues for the SAR model.
#'  String containing "L": Compute eigenvalues for the Leroux model.
#'  If a non-lattice model (spNNGP) is desired, set this to "None"; otherwise,
#'  match the lattice model desired.
#'
#' @return A list comprised of the following.
#' @returns coords_list: A list with a matrix of coordinates for each sample.
#' @returns covariates_list: A list with a dataframe of metadata for each sample.
#' @returns library_size_list: A list with a vector of library sizes for each sample.
#'  Total counts for each location if there are more than one gene/measurement;
#'  otherwise, set to 1.
#' @returns X_list: A list with a matrix of design matrices for each sample.
#' @returns W_list: A list with measurement adjacency matrix for each sample.
#'  If computed.
#' @returns D_list: A list with measurement degree matrix for each sample.
#'  If computed.
#' @returns eig_CS_list: A list with CAR/SAR model eigenvalues for each sample.
#'  If computed.
#' @returns eig_L_list: A list with Leroux model eigenvalues for each sample.
#'  If computed.
#'
#' @import Matrix
#' @importFrom Matrix t
#' @importFrom Matrix Diagonal
#' @importFrom Matrix rankMatrix
#' @importFrom Matrix rowSums
#' @importFrom Matrix solve
#' @importFrom Matrix sparseMatrix
#' @importFrom dplyr bind_rows
#' @importFrom spatstat.geom nnwhich
#' @importFrom spatstat.geom nndist
#' @importFrom stats model.matrix
#'
#' @export
prepData <- function(count_matrix,
                     meta_data,
                     sample_col,
                     design_formula = NULL,
                     design_mat = NULL,
                     coord_data = NULL,
                     D_THRESH = NULL,
                     k_search = 20,
                     adj_mat = NULL,
                     # C: CAR, S: SAR, L: Leroux, "NONE": None
                     compute_eigs = "CSL") {
  # Process arguments
  if (is.null(design_formula) && is.null(design_mat)) {
    stop("Must supply one of design formula or design matrix.")
  }
  if (!is.null(design_formula) && !is.null(design_mat)) {
    stop("Cannot supply both design formula and design matrix: pick one.")
  }
  if (is.null(compute_eigs)) {
    compute_eigs <- "NONE"
  }

  # Check dimensions and names
  stopifnot(ncol(count_matrix) == nrow(meta_data))
  stopifnot(identical(colnames(count_matrix), rownames(meta_data)))
  if (!is.null(coord_data)) {
    stopifnot(ncol(count_matrix) == nrow(coord_data))
    stopifnot(identical(rownames(coord_data), rownames(meta_data)))
  }
  if (!is.null(design_mat)) {
    stopifnot(ncol(count_matrix) == nrow(design_mat))
    stopifnot(identical(rownames(design_mat), rownames(meta_data)))
  }
  if (!is.null(adj_mat)) {
    cat("Using supplied adjacency matrix.", "\n")
    stopifnot(ncol(adj_mat) == nrow(adj_mat)) # Square
    stopifnot(ncol(adj_mat) == nrow(meta_data)) # Matching dimensions
    stopifnot(identical(rownames(adj_mat), colnames(adj_mat))) # Symmetry
    stopifnot(identical(rownames(adj_mat), rownames(meta_data))) # Names
  }

  # Sample names
  sample_names <- unique(meta_data[, sample_col])
  # Gene names
  gene_names <- rownames(count_matrix)

  # Lists to create for output
  # List of coordinate matrices
  coords_list <- list()
  # List of covariates data frames
  covariates_list <- list()
  # List of count matrices and library sizes
  counts_list <- list()
  library_size_list <- list()
  # List of design matrices
  X_list <- list()
  # List of adjacency matrices
  W_list <- list()
  D_list <- list()
  # Eigenvalues for CAR and SAR
  eig_CS_list <- list()
  # Eigenvalues for Leroux
  eig_L_list <- list()

  # Preprocess design_mat, if supplied
  if (!is.null(design_mat)) {
    cs_tmp <- colSums(design_mat)
    if (0 < sum(cs_tmp == 0)) {
      warning("There are all-zero columns in the design matrix.")
      cat("Which columns are all zero:", colnames(design_mat)[cs_tmp == 0], "\n")
      warning("Dropping columns that are all zero")
      design_mat <- design_mat[, cs_tmp != 0, drop = FALSE]
    }

    rank_tmp <- Matrix::rankMatrix(design_mat)[1]
    if (rank_tmp != ncol(design_mat)) {
      warning("The design matrix is not full rank.")
      cat("Design matrix is not full rank: ",
          rank_tmp,
          " < ",
          ncol(design_mat))
    }
  }

  # Extract the coordinates matrix for each sample
  # Extract the covariates to use for each sample
  # Extract the counts and library sizes for each sample
  for (samp in sample_names) {
    # Indices for cells in sample
    idx_local <- which(meta_data[, sample_col] == samp)
    cell_names_local <- rownames(meta_data)[idx_local]

    # Coordinates
    if (!is.null(coord_data)) {
      coords_list[[1 + length(coords_list)]] <- coord_data[idx_local, ]
      rownames(coords_list[[length(coords_list)]]) <- cell_names_local
    }

    # Covariates
    covariates_list[[1 + length(covariates_list)]] <- meta_data[idx_local, ]

    # Counts and library size
    counts_list[[1 + length(counts_list)]] <- count_matrix[, idx_local, drop =
                                                             FALSE]
    # Only use a library size that's the sum if we have more than one gene supplied
    if (1 < nrow(count_matrix)) {
      library_size_list[[1 + length(library_size_list)]] <- colSums(count_matrix[, idx_local, drop =
                                                                                   FALSE])
    } else {
      library_size_list[[1 + length(library_size_list)]] <- rep(1, length(cell_names_local))
    }


    # Design matrix
    if (!is.null(design_mat)) {
      X_list[[1 + length(X_list)]] <- design_mat[idx_local, , drop = FALSE]
      rownames(X_list[[length(X_list)]]) <- cell_names_local
    }

    # Set names
    names(library_size_list[[length(library_size_list)]]) <- cell_names_local
    rownames(counts_list[[length(counts_list)]]) <- gene_names
    colnames(counts_list[[length(counts_list)]]) <- cell_names_local
  }
  # Assign names to match samples
  names(coords_list) <- sample_names
  names(covariates_list) <- sample_names
  names(counts_list) <- sample_names
  names(library_size_list) <- sample_names

  ## Design matrix from metadata

  # Create design matrix if not supplied
  if (!is.null(design_formula)) {
    X_full <- dplyr::bind_rows(covariates_list)
    X_mat <- stats::model.matrix(design_formula, X_full)
    rm(X_full)

    # Check X_mat for zero columns and rank deficiency
    if (!is.null(X_mat)) {
      cs_tmp <- colSums(X_mat)
      if (0 < sum(cs_tmp == 0)) {
        warning("There are all-zero columns in the design matrix.")
        cat("Which columns are all zero:", colnames(X_mat)[cs_tmp == 0], "\n")
        warning("Dropping design matrix columns that are all zero")
        X_mat <- X_mat[, cs_tmp != 0, drop = FALSE]
      }

      rank_tmp <- Matrix::rankMatrix(X_mat)[1]
      if (rank_tmp != ncol(X_mat)) {
        warning("The design matrix is not full rank.")
        cat("Design matrix is not full rank: ",
            rank_tmp,
            " < ",
            ncol(X_mat))
      }
    }

    # Separate into lists
    for (samp in sample_names) {
      # Indices for cells in sample
      idx_local <- which(meta_data[, sample_col] == samp)
      cell_names_local <- rownames(meta_data)[idx_local]

      # Design matrix
      X_list[[1 + length(X_list)]] <- X_mat[idx_local, , drop = FALSE]
      rownames(X_list[[length(X_list)]]) <- cell_names_local
    }
  }
  names(X_list) <- sample_names

  # Adjacency matrix
  if (!is.null(adj_mat)) {
    cat("Subsetting provided adjacency matrix.", "\n")
    # Separate into lists
    for (samp in sample_names) {
      # Indices for cells in sample
      idx_local <- which(meta_data[, sample_col] == samp)
      cell_names_local <- rownames(meta_data)[idx_local]

      # Subset matrix
      W <- adj_mat[idx_local, idx_local]
      # Ensure symmetry (sometimes ties in distances lead to weird things)
      # W <- W + Matrix::t(W)
      # W <- 1 * (W != 0)
      # Check to make sure things are correct
      # stopifnot(max(abs(W - Matrix::t(W))) == 0)
      # stopifnot(max(W) == 1)
      # stopifnot(min(W) == 0)
      # stopifnot(0 == sum(abs(Matrix::diag(W)) != 0))

      # Degree matrix
      D <- Matrix::Diagonal(nrow(W), Matrix::rowSums(W))

      W_list[[1 + length(W_list)]] <- W
      D_list[[1 + length(D_list)]] <- D

      if (0 == min(Matrix::rowSums(W))) {
        warning(paste0("Not every cell has a neighbor, ", samp))
      }
      if (0 != max(abs(W - Matrix::t(W)))) {
        warning(paste0("Adjacency is not symmetric, ", samp))
      }
      if (1 != max(W)) {
        warning(paste0("Adjacency is not in {0, 1}, ", samp))
      }
      if (0 != min(W)) {
        warning(paste0("Adjacency is not in {0, 1}, ", samp))
      }
      if (0 != sum(abs(Matrix::diag(W)) != 0)) {
        warning(paste0("Adjacency is not diagonal-free, ", samp))
      }
    }
    names(W_list) <- sample_names
    names(D_list) <- sample_names
  } else if (!is.null(coord_data) &&
             !is.null(D_THRESH) && is.null(adj_mat)) {
    cat("Creating adjacency matrix.", "\n")

    ## Create adjacency matrices: If we have coordinates and a threshold
    for (idx in 1:length(coords_list)) {
      # Get indices of nearest neighbors for each cell
      nn_idx <- spatstat.geom::nnwhich(coords_list[[idx]], k = 1:k_search)
      # Get distances to nearest neighbors for each cell
      nn_dist <- spatstat.geom::nndist(coords_list[[idx]], k = 1:k_search)
      # Find number of neighbors for each cell--at least 1
      n_nb <- pmax(1, rowSums(nn_dist <= D_THRESH))

      # Extract indices of neighbors
      col_list <- list()
      row_list <- list()
      for (r_idx in 1:nrow(coords_list[[idx]])) {
        col_list[[r_idx]] <- nn_idx[r_idx, 1:n_nb[r_idx]]
        row_list[[r_idx]] <- rep(r_idx, n_nb[r_idx])
      }
      col_list <- Reduce(c, col_list)
      row_list <- Reduce(c, row_list)

      # Adjacency matrix of cell neighbors
      W <- Matrix::sparseMatrix(
        i = row_list,
        j = col_list,
        x = rep(1, length(row_list)),
        dims = c(nrow(coords_list[[idx]]), nrow(coords_list[[idx]]))
      )
      # Ensure symmetry (sometimes ties in distances lead to weird things)
      W <- W + Matrix::t(W)
      W <- 1 * (W != 0)
      # Check to make sure things are correct
      stopifnot(max(abs(W - Matrix::t(W))) == 0)
      stopifnot(max(W) == 1)
      stopifnot(min(W) == 0)
      stopifnot(0 == sum(abs(Matrix::diag(W)) != 0))

      # Degree matrix
      D <- Matrix::Diagonal(nrow(W), Matrix::rowSums(W))

      W_list[[idx]] <- W
      D_list[[idx]] <- D
    }
    names(W_list) <- sample_names
    names(D_list) <- sample_names
  } else {
    warning("Not creating adjacency matrices: missing either coordinates or threshold.")
  }

  ## Compute eigenvalues for the lattice models

  if (grepl("NONE", compute_eigs, ignore.case = TRUE)) {
    # Do nothing since we're told to not compute eigenvalues
    warning("Not computing eigenvalues: option selected.")
  } else if (0 == length(W_list)) {
    warning("Not computing eigenvalues: W_list is empty.")
  } else {
    for (idx in 1:length(W_list)) {
      cat(paste(Sys.time(), idx), "\n")
      if (grepl("C", compute_eigs, ignore.case = TRUE)
          || grepl("S", compute_eigs, ignore.case = TRUE)) {
        cat(
          "Starting CAR/SAR eigenvalue computation for area",
          idx,
          "at",
          paste(Sys.time(), "\n")
        )

        eig_CS_list[[idx]] <- Re(eigen(Matrix::solve(D_list[[idx]], W_list[[idx]]), FALSE, only.values = TRUE)$values)
      }
      if (grepl("L", compute_eigs, ignore.case = TRUE)) {
        cat(
          "Starting Leroux eigenvalue computation for area",
          idx,
          "at",
          paste(Sys.time(), "\n")
        )

        eig_L_list[[idx]] <- Re(eigen(D_list[[idx]] - W_list[[idx]], TRUE, only.values = TRUE)$values)
      }
      if (!(
        grepl("C", compute_eigs, ignore.case = TRUE)
        || grepl("S", compute_eigs, ignore.case = TRUE)
        || grepl("L", compute_eigs, ignore.case = TRUE)
      )) {
        warning("Not computing eigenvalues: Invalid options.")
      }
    }
    # Assign names
    if (grepl("C", compute_eigs, ignore.case = TRUE)
        || grepl("S", compute_eigs, ignore.case = TRUE)) {
      names(eig_CS_list) <- sample_names
    } else {
      eig_CS_list <- NULL
    }
    if (grepl("L", compute_eigs, ignore.case = TRUE)) {
      names(eig_L_list) <- sample_names
    } else {
      eig_L_list <- NULL
    }
  }

  ## Store and return
  return(structure(
    list(
      coords_list = coords_list,
      covariates_list = covariates_list,
      counts_list = counts_list,
      library_size_list = library_size_list,
      X_list = X_list,
      W_list = W_list,
      D_list = D_list,
      eig_CS_list = eig_CS_list,
      eig_L_list = eig_L_list
    ),
    class = "TESSERAData"
  ))
}

#' Prepare data for the TESSERA method: Choose a distance threshold for
#'  adjacency matrices.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param meta_data Dataframe with metadata/covariates.
#'  Ensure that rownames(meta_data) == rownames(coord_data).
#' @param sample_col String: Column name in meta_data identifying which rows correspond
#'  to which sample (how to break up the data into multiple samples).
#' @param coord_data Dataframe or matrix with coordinates for observations.
#'  Rows are samples.
#'  Ensure that rownames(coord_data) == rownames(meta_data).
#' @param k_search When forming adjacency matrices, the maximum number of neighbors.
#'
#' @return A list comprised of the following.
#' @returns nb_dist: Dataframe containing distances to the k-th nearest neighbor
#'  for each cell for each sample.
#'  First column is sample, all other columns are the distances to the k-th nearest neighbor.
#' @returns mean_nb_dist: Matrix containing mean distance (across cells) to
#'  k-th nearest neighbor for each sample.
#'  Rows are samples, columns are mean distances to the kth nearest neighbor.
#' @returns std_nb_dist: Matrix containing standard deviation of distance (across cells) to
#'  k-th nearest neighbor for each sample.
#'  Rows are samples, columns are standard deviation of distances to the kth nearest neighbor.
#' @returns distances_plot: ggplot2 plot showing mean Euclidean distance (across cells)
#'  for each sample to the k-th nearest neighbor.
#'
#' @import ggplot2
#' @importFrom dplyr bind_rows
#' @importFrom spatstat.geom nndist
#' @importFrom reshape2 melt
#' @importFrom Rfast colVars
#'
#' @export
visualizeNeighborDistances <- function(meta_data,
                                       sample_col,
                                       coord_data,
                                       k_search = 20) {
  ## List of spatial coordinates for each sample
  coords_list <- lapply(unique(meta_data[, sample_col]), function (x) {
    coord_data[meta_data[, sample_col] == x, ]
  })
  names(coords_list) <- unique(meta_data[, sample_col])

  ## k nearest neighbor distances for each cell
  nb_dist <- dplyr::bind_rows(lapply(1:length(coords_list), function (x) {
    tmp <- spatstat.geom::nndist(coords_list[[x]], k = 1:k_search)
    rownames(tmp) <- rownames(coords_list[[x]])
    tmp <- cbind(rep(names(coords_list)[x], nrow(tmp)), tmp)
    colnames(tmp)[1] = "sample"
    data.frame(tmp)
  }))

  ## Mean / Standard deviation of k-nearest neighbor distances
  # Each row is a sample, each column is a distance (compute metrics over cells)
  mean_nb_dist <- as.matrix(dplyr::bind_rows(lapply(coords_list, function (x) {
    colMeans(spatstat.geom::nndist(x, k = 1:k_search))
  })))
  rownames(mean_nb_dist) <- names(coords_list)

  std_nb_dist <- sqrt(t(as.matrix(dplyr::bind_cols(
    lapply(coords_list, function (x) {
      as.matrix(Rfast::colVars(spatstat.geom::nndist(x, k = 1:k_search)))
    })
  ))))
  rownames(std_nb_dist) <- names(coords_list)

  ## Make a plot to visualize mean Euclidean distance of cells to kth neighbor for each sample
  plot_data <- reshape2::melt(mean_nb_dist,
                              value.name = "distance",
                              varnames = c("index", "k_nghbr"))

  plot_data$k_nghbr <- gsub("dist.", "", plot_data$k_nghbr)
  plot_data$k_nghbr <- factor(plot_data$k_nghbr, levels = 1:k_search)
  p <- ggplot2::ggplot(data = plot_data, ggplot2::aes_string(x = "k_nghbr", y = "distance"))
  p <- p + ggplot2::geom_boxplot()
  p <- p + ggplot2::geom_jitter()
  p <- p + ggplot2::labs(x = "kth neighbor", y = "Mean Euclidean distance")
  p <- p + ggplot2::theme_bw()
  p <- p + ggplot2::theme(text = ggplot2::element_text(size = 20))

  ## Store and return
  return(
    list(
      nb_dist = nb_dist,
      mean_nb_dist = mean_nb_dist,
      std_nb_dist = std_nb_dist,
      distances_plot = p
    )
  )
}
