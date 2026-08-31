#' The TESSERAOutput Class
#'
#' The \code{TESSERAOutput} object is the primary model-fitting output container
#' for the \code{TESSERA} package, returned by \code{\link{TESSERA_lattice}} and
#' \code{\link{TESSERA_spNNGP}}. It stores estimated fixed effects, sample-specific
#' spatial random effects, spatial covariance parameters, optimization histories,
#' goodness-of-fit diagnostics, and negative Hessian matrices for downstream inference.
#'
#' Downstream functions—-specifically \code{\link{calc_Wald_statistics}} (for linear
#' contrast testing) and \code{\link{summarize_TESSERA}} (for sample-level performance
#' summaries)—-ingest \code{TESSERAOutput} objects directly.
#' 
#' @section Elements:
#' A \code{TESSERAOutput} object is a structured list containing the following components:
#' \itemize{
#'   \item \strong{beta_hat}: Numeric vector of estimated shared fixed effect coefficients \eqn{\hat{\boldsymbol{\beta}}}.
#'   \item \strong{gamma_hat}: Numeric vector of estimated spatial correlation/dependence parameters \eqn{\hat{\gamma}_i} for each sample (\emph{lattice models only}).
#'   \item \strong{tau2_hat}: Numeric vector of estimated spatial variance scale parameters \eqn{\hat{\tau}_i^2} for each sample (\emph{lattice models only}).
#'   \item \strong{cov_param_hat}: Numeric matrix of estimated spatial covariance parameters (Nugget, Sill, Range, Smoothness) across samples (\emph{spNNGP models only}).
#'   \item \strong{phi_hat}: Numeric vector of estimated spatial random effects \eqn{\hat{\boldsymbol{\phi}} = \hat{\boldsymbol{\eta}} - \boldsymbol{X}\hat{\boldsymbol{\beta}}} across all observations.
#'   \item \strong{eta_hat}: Numeric vector of latent linear predictors \eqn{\hat{\boldsymbol{\eta}}} (E-step conditional expectation means).
#'   \item \strong{theta_hat}: Numeric vector of estimated Poisson relative rate parameters \eqn{\exp(\boldsymbol{X}\hat{\boldsymbol{\beta}} + \hat{\boldsymbol{\phi}})}.
#'   \item \strong{predictions}: Numeric vector of predicted counts \eqn{\hat{z} = \hat{\theta} \times \text{offset}} across all observations.
#'   \item \strong{residuals}: Numeric vector of raw residuals (\eqn{z - \hat{z}}) across all observations.
#'   \item \strong{beta_tracker}: Numeric matrix recording the trajectory of \eqn{\boldsymbol{\beta}} estimates across ECM iterations (parameters \eqn{\times} iterations).
#'   \item \strong{gamma_tracker} / \strong{cov_param_tracker}: History of spatial parameter estimates across iterations.
#'   \item \strong{tau2_tracker}: Numeric matrix recording the trajectory of \eqn{\tau^2} estimates across iterations (\emph{lattice models only}).
#'   \item \strong{R2_tracker}: Numeric matrix of sample-wise \eqn{R^2} correlations across iterations.
#'   \item \strong{MSE_tracker}: Numeric matrix of sample-wise mean squared errors across iterations.
#'   \item \strong{data_log_like_tracker}: Numeric matrix of observed Poisson data log-likelihoods across iterations.
#'   \item \strong{expected_log_like_tracker}: Numeric matrix of expected complete-data log-likelihoods (\eqn{Q}-function) across iterations.
#'   \item \strong{resid_moran}: Array/matrix recording residual spatial autocorrelation (Moran's \eqn{I}, expected value, and \eqn{p}-value) across iterations.
#'   \item \strong{resid_moran_nb}: Matrix recording neighborhood-based residual Moran's \eqn{I} across iterations (\emph{lattice models only}).
#'   \item \strong{beta_neghessian}: Negative Hessian matrix evaluated at \eqn{\hat{\boldsymbol{\beta}}}, representing the observed Fisher information used to compute standard errors and Wald statistics.
#'   \item \strong{tau2_neghessian}: Vector of negative Hessians for \eqn{\tau^2} per sample (\emph{lattice models only}).
#'   \item \strong{gamma_neghessian}: Vector of negative Hessians for \eqn{\gamma} per sample (\emph{lattice models only}).
#'   \item \strong{performanceSummary}: Data frame summarizing sample-wise model performance, goodness-of-fit (normalized MSE, Moran's \eqn{I}), and parameter estimates.
#'   \item \strong{start_idx_list}: Named integer vector indicating the starting index of each sample in the concatenated observation vectors.
#'   \item \strong{time}: \code{difftime} object recording the total elapsed optimization time.
#'   \item \strong{run_settings}: List of runtime configuration options, convergence thresholds, and initialization flags used during fitting.
#' }
#'
#' @seealso \code{\link{TESSERA_lattice}}, \code{\link{TESSERA_spNNGP}}, \code{\link{calc_Wald_statistics}}, \code{\link{summarize_TESSERA}}
#'
#' @name TESSERAOutput
#' @aliases TESSERAOutput-class
NULL
