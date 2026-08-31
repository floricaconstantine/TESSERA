#' TESSERA: Tool for Estimating Spatial and Sample-level Effects via Regression Analysis
#'
#' The \code{TESSERA} package implements multi-sample spatial generalized linear
#' mixed models (GLMM) for spatial count data across disconnected samples. While
#' primarily developed for differential expression in spatial transcriptomics, it
#' provides a general statistical framework for estimating a shared set of fixed
#' effects across multi-sample designs while accommodating independent,
#' sample-specific spatial covariance structures and distinct coordinate systems.
#' 
#'
#' Key features of \code{TESSERA} include:
#' \itemize{
#'   \item \strong{Multi-Sample GLMMs:} Models overdispersed Poisson count data with sample-specific spatial random effects and a shared set of fixed effects across all samples.
#'   \item \strong{Flexible Covariance Structures:} Supports both lattice models (Leroux, CAR, SAR) and sparse nearest-neighbor Gaussian processes (spNNGP with Matérn, Exponential, Gaussian, or Spherical kernels).
#'   \item \strong{Efficient Parameter Estimation:} Employs an Expectation-Conditional-Maximization (ECM) algorithm with pre-computed spatial matrix operations for scalable gene-by-gene fitting.
#'   \item \strong{Flexible Hypothesis Testing:} Computes effect size estimates and Wald statistics for user-specified linear contrasts to test diverse scientific hypotheses (such as subgroup-specific treatment effects or cell-type comparisons).
#'   \item \strong{Empirical Null Calibration:} Estimates a scaled non-central \eqn{\chi_1^2} empirical null distribution to correct for finite-sample GLMM bias and control the False Discovery Rate (FDR).
#' }
#'
#' @author Florica Constantine \email{florica@berkeley.edu}
#'
#' @references
#' Constantine, F. et al. (2026). Unlocking Multi-Sample Differential Expression
#' for Spatial Transcriptomics Data with TESSERA. \emph{bioRxiv}.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{prep_data}}: Assembles input data matrices or \code{SpatialExperiment} objects and pre-computes spatial neighborhood structures and eigenvalues.
#'   \item \code{\link{plot_neighbor_distances}}: Diagnostic tool to visualize nearest-neighbor distances and determine adjacency thresholds for lattice models.
#'   \item \code{\link{TESSERA_lattice}}: Fits an overdispersed Poisson spatial GLMM for a single gene using lattice covariance structures (Leroux, CAR, or SAR).
#'   \item \code{\link{TESSERA_spNNGP}}: Fits an overdispersed Poisson spatial GLMM for a single gene using spNNGP covariance kernels (Matérn, Exponential, Gaussian, or Spherical kernels).
#'   \item \code{\link{calc_Wald_statistics}}: Computes Wald test statistics, contrast estimates, and standard errors for user-specified contrast matrices.
#'   \item \code{\link{select_Wald_threshold}}: Sweeps over candidate thresholds to fit the scaled non-central \eqn{\chi_1^2} empirical null distribution.
#'   \item \code{\link{calc_scaled_noncentral_chi2_pvalues}}: Calculates empirical \eqn{p}-values from Wald statistics using the fitted empirical null parameters.
#' }
#'
#' @keywords internal
"_PACKAGE"
