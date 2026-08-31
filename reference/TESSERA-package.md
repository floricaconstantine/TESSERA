# TESSERA: Tool for Estimating Spatial and Sample-level Effects via Regression Analysis

The `TESSERA` package implements multi-sample spatial generalized linear
mixed models (GLMM) for spatial count data across disconnected samples.
While primarily developed for differential expression in spatial
transcriptomics, it provides a general statistical framework for
estimating a shared set of fixed effects across multi-sample designs
while accommodating independent, sample-specific spatial covariance
structures and distinct coordinate systems.

## Details

Key features of `TESSERA` include:

- **Multi-Sample GLMMs:** Models overdispersed Poisson count data with
  sample-specific spatial random effects and a shared set of fixed
  effects across all samples.

- **Flexible Covariance Structures:** Supports both lattice models
  (Leroux, CAR, SAR) and sparse nearest-neighbor Gaussian processes
  (spNNGP with Matérn, Exponential, Gaussian, or Spherical kernels).

- **Efficient Parameter Estimation:** Employs an
  Expectation-Conditional-Maximization (ECM) algorithm with pre-computed
  spatial matrix operations for scalable gene-by-gene fitting.

- **Flexible Hypothesis Testing:** Computes effect size estimates and
  Wald statistics for user-specified linear contrasts to test diverse
  scientific hypotheses (such as subgroup-specific treatment effects or
  cell-type comparisons).

- **Empirical Null Calibration:** Estimates a scaled non-central
  \\\chi_1^2\\ empirical null distribution to correct for finite-sample
  GLMM bias and control the False Discovery Rate (FDR).

## Main Functions

- [`prep_data`](https://floricaconstantine.github.io/TESSERA/reference/prep_data.md):
  Assembles input data matrices or `SpatialExperiment` objects and
  pre-computes spatial neighborhood structures and eigenvalues.

- [`plot_neighbor_distances`](https://floricaconstantine.github.io/TESSERA/reference/plot_neighbor_distances.md):
  Diagnostic tool to visualize nearest-neighbor distances and determine
  adjacency thresholds for lattice models.

- [`TESSERA_lattice`](https://floricaconstantine.github.io/TESSERA/reference/TESSERA_lattice.md):
  Fits an overdispersed Poisson spatial GLMM for a single gene using
  lattice covariance structures (Leroux, CAR, or SAR).

- [`TESSERA_spNNGP`](https://floricaconstantine.github.io/TESSERA/reference/TESSERA_spNNGP.md):
  Fits an overdispersed Poisson spatial GLMM for a single gene using
  spNNGP covariance kernels (Matérn, Exponential, Gaussian, or Spherical
  kernels).

- [`calc_Wald_statistics`](https://floricaconstantine.github.io/TESSERA/reference/calc_Wald_statistics.md):
  Computes Wald test statistics, contrast estimates, and standard errors
  for user-specified contrast matrices.

- [`select_Wald_threshold`](https://floricaconstantine.github.io/TESSERA/reference/select_Wald_threshold.md):
  Sweeps over candidate thresholds to fit the scaled non-central
  \\\chi_1^2\\ empirical null distribution.

- [`calc_scaled_noncentral_chi2_pvalues`](https://floricaconstantine.github.io/TESSERA/reference/calc_scaled_noncentral_chi2_pvalues.md):
  Calculates empirical \\p\\-values from Wald statistics using the
  fitted empirical null parameters.

## References

Constantine, F. et al. (2026). Unlocking Multi-Sample Differential
Expression for Spatial Transcriptomics Data with TESSERA. *bioRxiv*.

## See also

Useful links:

- <https://floricaconstantine.github.io/TESSERA/>

## Author

Florica Constantine <florica@berkeley.edu>
