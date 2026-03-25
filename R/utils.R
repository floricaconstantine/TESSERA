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
#' ex_data <- generate_data_one_area(1000, 0.03, "Leroux", beta_true, gamma_true, tau2_true, "rand_bern")
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
#' Taken from: https://stackoverflow.com/questions/5560218/computing-sparse-pairwise-distance-matrix-in-r
#'
#' @param coords (x, y) for points as rows.
#' @param k Number of nearest neighbors.
#'
#' @returns 2 k x n points matrix.
#'  First k rows are distances, last k are indices of nearest neighbors.
#'  Only Upper Triangular part is formed.
sparseDist <- function(coords, k) {
  # Transpose coordinates so that columns are samples
  coords <- t(coords)
  # Number of samples
  n <- ncol(coords)

  d <- vapply(seq_len(n - 1L), function(i) {
    # Squared L2 distances to all points with higher indices
    d <- colSums((coords[, seq(i + 1L, n), drop = FALSE] - coords[, i])^2)
    # Sort distances and subset to top k
    o <- sort.list(d, na.last = NA, method = "quick")[seq_len(k)]
    # Return vector of distances (L2) and indices
    c(sqrt(d[o]), o + i)
  }, numeric(2 * k))
  dimnames(d) <- list(c(paste('d', seq_len(k), sep = ''), paste('i', seq_len(k), sep =
                                                                  '')), colnames(coords)[-n])
  return(d)
}

#' Get k nearest neighbors/distances given coordinates.
#' Modified from: https://stackoverflow.com/questions/5560218/computing-sparse-pairwise-distance-matrix-in-r
#'
#' @param coords (x, y) for points as rows.
#' @param k Number of nearest neighbors.
#'
#' @returns 2 k x n points matrix.
#'  First k rows are distances, last k are indices of nearest neighbors.
#'  Only Lower Triangular part is formed.
sparseDist_LT <- function(coords, k) {
  # Transpose coordinates so that columns are samples
  coords <- t(coords)
  # Number of samples
  n <- ncol(coords)

  d <- vapply(seq(2, n), function(i) {
    # Squared L2 distances to all points with lower indices
    d <- colSums((coords[, seq(1, i - 1), drop = FALSE] - coords[, i])^2)
    # Sort distances and subset to top k
    o <- sort.list(d, na.last = NA, method = "quick")[seq_len(k)]
    # Return vector of distances (L2) and indices
    c(sqrt(d[o]), o)
  }, numeric(2 * k))
  dimnames(d) <- list(c(paste('d', seq_len(k), sep = ''), paste('i', seq_len(k), sep =
                                                                  '')), colnames(coords)[-n])
  return(d)
}
