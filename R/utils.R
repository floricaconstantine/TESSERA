## General utilties that are not model or algorithm specific.
# Dependencies in file: N/A.

#' Compute Moran's I using an adjacency matrix.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param y A vector of measurements.
#' @param W Neighbor/adjacency matrix (symmetric, binary).
#'
#' @returns Moran's I (spatial autocorrelation).
#' @export
#'
#' @examples
#' set.seed(2026)
#' tau2_true <- 1.0
#' gamma_true <- 0.5
#' beta_true <- c(1, 0, -1)
#' ex_data <- generate_data_one_area(1000, 0.03, "Leroux", beta_true,
#'   gamma_true, tau2_true, "rand_bern")
#' moran_I_nb(ex_data$z, ex_data$W)
moran_I_nb <- function(y, W) {
  # Subtract mean
  ym <- y - mean(y)
  
  # Numerator: Quadratic form
  MI <- (t(ym) %*% W) %*% ym
  # Denominator: Sum of squares
  MI <- MI / sum(ym * ym)
  # Scaling
  MI <- MI * (length(ym) / sum(W))
  # Ensure is a number not a 1x1 matrix
  MI <- as.numeric(MI)
  
  return(MI)
}

#' Compute the log likelihood of a set of Poisson variables.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Discards invalid values in computation.
#' @note Discards factorial terms for numerical stability, etc.
#'
#' @param z Observed counts.
#' @param theta_hat Estimated Poisson parameters.
#'
#' @returns log likelihood of data.
poisson_loglike <- function(z, theta_hat) {
  # Data log likelihood
  # tmp_ll <- dpois(z,
  #                 pmax(0.0, theta_hat),
  #                 log = TRUE)
  tmp_ll <- pmax(0.0, z) * log(pmax(0.0, theta_hat)) - pmax(0.0, theta_hat)
  tmp_ll[is.infinite(tmp_ll)] <- 0.0
  tmp_ll[is.nan(tmp_ll)] <- 0.0
  
  return(sum(tmp_ll))
}

#' Utility function to instantiate variables from elements in a named list.
#'  E.g., call this function on list(a=1, b=2) would result in variables a and b
#'  with values 1 and 2, respectively, in the calling frame/environment.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param lst List with named fields.
#' @param target_environ Environment to instantiate variables in.
#'
#' @note This function is not currently used elsewhere in the codebase.
#' @note This function has identical functionality and essentially identical code
#'  to extract.named in mvbutils.
#'
#' @returns Nothing; the global environment is modified.
variables_from_list <- function(lst, target_environ = parent.frame()) {
  # Get names of fields in list lst
  list_fields <- names(lst)
  # Loop over fields with names
  for (fld in list_fields[nchar(list_fields) > 0]) {
    assign(fld, lst[[fld]], envir = target_environ)
  }
}

#' Get k nearest neighbors/distances given coordinates.
#'
#' For each point, we only look at points with a higher index.
#' E.g., for point 10, we only look at points 11, 12, ...
#' We do not look at the last point (since there are no higher-index points).
#' So the output has one fewer column than the number of points.
#' If there are fewer than k points with higher index, we fill the remaining
#' entries with NA.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param coords (x, y) for points as rows.
#' @param k Number of nearest neighbors.
#'
#' @returns 2 k x n points matrix.
#'  Top k rows: Euclidean distances to nearest neighbors.
#'  Bottom k rows: Indices of nearest neighbors.
#'  Each column corresponds to a point.
#' 
#' @export
#' 
#' @examples
#' sparse_dist(matrix(1:4, ncol=2), 1)
sparse_dist <- function(coords, k) {
  coords_mat <- if (is.matrix(coords))
    coords
  else
    as.matrix(coords)
  storage.mode(coords_mat) <- "double"
  return(sparse_dist_cpp(coords_mat, as.integer(k)))
}

#' Get k nearest neighbors/distances given coordinates.
#'
#' For each point, we only look at points with a lower index.
#' E.g., for point 10, we only look at points 1, 2, ..., 9.
#' We do not look at the first point (since there are no lower-index points).
#' So the output has one fewer column than the number of points.
#' If there are fewer than k points with lower index, we fill the remaining
#' entries with NA.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @note Applies to a single area.
#'
#' @param coords (x, y) for points as rows.
#' @param k Number of nearest neighbors.
#'
#' @returns 2 k x n points matrix.
#'  Top k rows: Euclidean distances to nearest neighbors.
#'  Bottom k rows: Indices of nearest neighbors.
#'  Each column corresponds to a point.
#'  
#' @export
#' 
#' @examples
#' sparse_dist_LT(matrix(1:4, ncol=2), 1)
sparse_dist_LT <- function(coords, k) {
  coords_mat <- if (is.matrix(coords))
    coords
  else
    as.matrix(coords)
  storage.mode(coords_mat) <- "double"
  return(sparse_dist_LT_cpp(coords_mat, as.integer(k)))
}
