# Package index

## All functions

- [`BRISC_wrapper()`](https://floricaconstantine.github.io/TESSERA/reference/BRISC_wrapper.md)
  : Wrapper to fit sparse NN GP model via BRISC.
- [`E_step_Vhat()`](https://floricaconstantine.github.io/TESSERA/reference/E_step_Vhat.md)
  : Compute the Covariance Matrix of the random effects eta. Part of the
  E-Step in the EM algorithm. Follows Clayton and Kaldor (1987) for the
  Poisson-Gaussian model.
- [`E_step_etahat()`](https://floricaconstantine.github.io/TESSERA/reference/E_step_etahat.md)
  : Compute the Expectation of the random effects eta. Part of the
  E-Step in the EM algorithm. Follows Clayton and Kaldor (1987) for the
  Poisson-Gaussian model.
- [`E_step_predict()`](https://floricaconstantine.github.io/TESSERA/reference/E_step_predict.md)
  : Compute the best fit to the data. Associated with the E-Step in the
  EM algorithm, but not needed to run. I.e., this is a utility function.
  Follows Clayton and Kaldor (1987) for the Poisson-Gaussian model.
- [`E_step_thetahat()`](https://floricaconstantine.github.io/TESSERA/reference/E_step_thetahat.md)
  : Compute the estimated poisson parameters theta. Note that eta = log
  theta. Associated with the E-Step in the EM algorithm, but not needed
  to run. I.e., this is a utility function. Follows Clayton and
  Kaldor (1987) for the Poisson-Gaussian model.
- [`M_step_BRISC()`](https://floricaconstantine.github.io/TESSERA/reference/M_step_BRISC.md)
  : Holding beta constant, use BRISC to find the spatial parameters of
  eta - X beta. Part of the M-Step in the EM algorithm. ONLY APPLIES TO
  THE GAUSSIAN PROCESS MODEL.
- [`M_step_beta()`](https://floricaconstantine.github.io/TESSERA/reference/M_step_beta.md)
  : Maximize the expected likelihood in beta, holding other variables
  constant. Part of the M-Step in the EM algorithm.
- [`M_step_gamma_CAR()`](https://floricaconstantine.github.io/TESSERA/reference/M_step_gamma_CAR.md)
  : Maximize the expected likelihood in gamma, holding other variables
  constant. Part of the M-Step in the EM algorithm. ONLY APPLIES TO THE
  CAR MODEL.
- [`M_step_gamma_Leroux()`](https://floricaconstantine.github.io/TESSERA/reference/M_step_gamma_Leroux.md)
  : Maximize the expected likelihood in gamma, holding other variables
  constant. Part of the M-Step in the EM algorithm. ONLY APPLIES TO THE
  LEROUX MODEL.
- [`M_step_gamma_SAR()`](https://floricaconstantine.github.io/TESSERA/reference/M_step_gamma_SAR.md)
  : Maximize the expected likelihood in gamma, holding other variables
  constant. Part of the M-Step in the EM algorithm. ONLY APPLIES TO THE
  SAR MODEL.
- [`M_step_tau2()`](https://floricaconstantine.github.io/TESSERA/reference/M_step_tau2.md)
  : Maximize the expected likelihood in tau^2, holding other variables
  constant. Part of the M-Step in the EM algorithm.
- [`M_step_variogram()`](https://floricaconstantine.github.io/TESSERA/reference/M_step_variogram.md)
  : Holding beta constant, fit a variogram model to eta - X beta. Part
  of the M-Step in the EM algorithm. Instead of a traditional MLE, which
  has cubic time complexity, we fit a variogram to estimate the kernel
  parameters, which is quadratic time. ONLY APPLIES TO THE GAUSSIAN
  PROCESS MODEL.
- [`Q_matrix_CAR()`](https://floricaconstantine.github.io/TESSERA/reference/Q_matrix_CAR.md)
  : Compute unscaled precision matrix for CAR model.
- [`Q_matrix_Leroux()`](https://floricaconstantine.github.io/TESSERA/reference/Q_matrix_Leroux.md)
  : Compute unscaled precision matrix for Leroux model.
- [`Q_matrix_SAR()`](https://floricaconstantine.github.io/TESSERA/reference/Q_matrix_SAR.md)
  : Compute unscaled precision matrix for SAR model.
- [`TESSERAData`](https://floricaconstantine.github.io/TESSERA/reference/TESSERAData.md)
  [`TESSERAData-class`](https://floricaconstantine.github.io/TESSERA/reference/TESSERAData.md)
  : The TESSERAData Class
- [`TESSERAOutput`](https://floricaconstantine.github.io/TESSERA/reference/TESSERAOutput.md)
  [`TESSERAOutput-class`](https://floricaconstantine.github.io/TESSERA/reference/TESSERAOutput.md)
  : The TESSERAOutput Class
- [`TESSERA_lattice()`](https://floricaconstantine.github.io/TESSERA/reference/TESSERA_lattice.md)
  : Fit Multi-Sample Poisson Spatial GLMM via ECM Algorithm
- [`TESSERA_spNNGP()`](https://floricaconstantine.github.io/TESSERA/reference/TESSERA_spNNGP.md)
  : Fit Multi-Sample Poisson Spatial GLMM via spNNGP
- [`build_Vinv_cpp()`](https://floricaconstantine.github.io/TESSERA/reference/build_Vinv_cpp.md)
  : Fast C++ Construction of E-Step Precision Matrix Vinv = Q / tau2 +
  diag(diag_add).
- [`build_nb_dist_cpp()`](https://floricaconstantine.github.io/TESSERA/reference/build_nb_dist_cpp.md)
  : Precompute Static Euclidean Distance Matrices for Neighbors in C++.
- [`calc_Wald_pvalue_from_fdrtool()`](https://floricaconstantine.github.io/TESSERA/reference/calc_Wald_pvalue_from_fdrtool.md)
  : Given Wald statistics, compute p-values.
- [`calc_Wald_pvalue_from_threshold()`](https://floricaconstantine.github.io/TESSERA/reference/calc_Wald_pvalue_from_threshold.md)
  : Given Wald statistics, compute p-values.
- [`calc_Wald_statistics()`](https://floricaconstantine.github.io/TESSERA/reference/calc_Wald_statistics.md)
  : Given the output of the TESSERA algorithms and a contrast matrix,
  compute Wald T-statistics.
- [`calc_Wald_statistics_old()`](https://floricaconstantine.github.io/TESSERA/reference/calc_Wald_statistics_old.md)
  : Given the output of the TESSERA algorithms and a contrast matrix,
  compute Wald T-statistics.
- [`calc_moran()`](https://floricaconstantine.github.io/TESSERA/reference/calc_moran.md)
  : Fast computation of Moran's I.
- [`calc_moran_batch_cpp()`](https://floricaconstantine.github.io/TESSERA/reference/calc_moran_batch_cpp.md)
  : Batch Computation of Moran's I Across Multiple Vectors on Shared
  Coordinates.
- [`calc_moran_precomp()`](https://floricaconstantine.github.io/TESSERA/reference/calc_moran_precomp.md)
  : Fast computation of Moran's I using precomputed coordinate
  constants.
- [`calc_scaled_noncentral_chi2_pvalues()`](https://floricaconstantine.github.io/TESSERA/reference/calc_scaled_noncentral_chi2_pvalues.md)
  : Given Wald statistics, compute p-values.
- [`check_inputs_TESSERA()`](https://floricaconstantine.github.io/TESSERA/reference/check_inputs_TESSERA.md)
  : Check inputs for the TESSERA method. Make sure that the input object
  has everything needed to run without error.
- [`check_inputs_TESSERA_spNNGP()`](https://floricaconstantine.github.io/TESSERA/reference/check_inputs_TESSERA_spNNGP.md)
  : Check inputs for the TESSERA_spNNGP method. Make sure that the input
  object has everything needed to run without error.
- [`distanceCalculate()`](https://floricaconstantine.github.io/TESSERA/reference/distanceCalculate.md)
  : Helper function for Moran's I.
- [`eval_Q_lattice_cpp()`](https://floricaconstantine.github.io/TESSERA/reference/eval_Q_lattice_cpp.md)
  : Fast C++ Evaluation of Lattice Precision Matrix Q from Pre-aligned
  CSC Vectors.
- [`example_TESSERA_out_Leroux.rds`](https://floricaconstantine.github.io/TESSERA/reference/example_TESSERA_out_Leroux.rds.md)
  : Fitted TESSERA Leroux Model Object
- [`example_prepData.rds`](https://floricaconstantine.github.io/TESSERA/reference/example_prepData.rds.md)
  : Pre-processed TESSERA Data Object
- [`example_raw_data.rds`](https://floricaconstantine.github.io/TESSERA/reference/example_raw_data.rds.md)
  : Raw Synthetic Spatial Data
- [`expected_loglike()`](https://floricaconstantine.github.io/TESSERA/reference/expected_loglike.md)
  : Compute the expected log likelihood.
- [`fast_trace_quad_cpp()`](https://floricaconstantine.github.io/TESSERA/reference/fast_trace_quad_cpp.md)
  : Fast Trace of Product of Two Sparse CSC Matrices Tr(Q \* Vhat).
- [`fit_empirical_null_fdrtool()`](https://floricaconstantine.github.io/TESSERA/reference/fit_empirical_null_fdrtool.md)
  : Wrapper for Empirical null estimation and p-value computation
- [`fit_scaled_noncentral_chi2()`](https://floricaconstantine.github.io/TESSERA/reference/fit_scaled_noncentral_chi2.md)
  : Fit a scaled non-central chi^2_1 distribution using BOBYQA
- [`fit_scaled_noncentral_chi2_old()`](https://floricaconstantine.github.io/TESSERA/reference/fit_scaled_noncentral_chi2_old.md)
  : Fit a scaled non-central chi^2_1 distribution using BOBYQA
- [`generate_data_one_area()`](https://floricaconstantine.github.io/TESSERA/reference/generate_data_one_area.md)
  : Simulate spatial count data for a single area
- [`generate_data_one_area_spNNGP()`](https://floricaconstantine.github.io/TESSERA/reference/generate_data_one_area_spNNGP.md)
  : Simulate spatial count data via spNNGP
- [`invert_precision_matrix()`](https://floricaconstantine.github.io/TESSERA/reference/invert_precision_matrix.md)
  : Helper function for Wald stat covariances.
- [`kernel_exp()`](https://floricaconstantine.github.io/TESSERA/reference/kernel_exp.md)
  : Exponential Kernel Function.
- [`kernel_gauss()`](https://floricaconstantine.github.io/TESSERA/reference/kernel_gauss.md)
  : Gaussian (Squared Exponential) Kernel Function.
- [`kernel_matern()`](https://floricaconstantine.github.io/TESSERA/reference/kernel_matern.md)
  : Matern Kernel Function.
- [`kernel_sph()`](https://floricaconstantine.github.io/TESSERA/reference/kernel_sph.md)
  : Spherical Kernel Function.
- [`moran_I_nb()`](https://floricaconstantine.github.io/TESSERA/reference/moran_I_nb.md)
  : Compute Moran's I using an adjacency matrix.
- [`neg_hessian_beta()`](https://floricaconstantine.github.io/TESSERA/reference/neg_hessian_beta.md)
  : Compute the Negative Hessian for beta.
- [`neg_hessian_gamma_CAR()`](https://floricaconstantine.github.io/TESSERA/reference/neg_hessian_gamma_CAR.md)
  : Compute the Negative Hessian of gamma. ONLY APPLIES TO THE CAR
  MODEL.
- [`neg_hessian_gamma_Leroux()`](https://floricaconstantine.github.io/TESSERA/reference/neg_hessian_gamma_Leroux.md)
  : Compute the Negative Hessian of gamma. ONLY APPLIES TO THE Leroux
  MODEL.
- [`neg_hessian_gamma_SAR()`](https://floricaconstantine.github.io/TESSERA/reference/neg_hessian_gamma_SAR.md)
  : Compute the Negative Hessian of gamma. ONLY APPLIES TO THE SAR
  MODEL.
- [`neg_hessian_tau2()`](https://floricaconstantine.github.io/TESSERA/reference/neg_hessian_tau2.md)
  : Compute the Negative Hessian for tau^2.
- [`nngp_beta_terms_cpp()`](https://floricaconstantine.github.io/TESSERA/reference/nngp_beta_terms_cpp.md)
  : Compute X^T Q eta and X^T Q X Directly from NNGP Factors (A, Dinv)
  Without Forming Q.
- [`nngp_prec_mat()`](https://floricaconstantine.github.io/TESSERA/reference/nngp_prec_mat.md)
  : Form a sparse Nearest-Neighbor Gaussian Process Precision Matrix.
- [`nngp_prec_triplets_cpp()`](https://floricaconstantine.github.io/TESSERA/reference/nngp_prec_triplets_cpp.md)
  : Compute Sparse NNGP Precision Factors (A triplets and Dinv) in C++.
- [`normalize()`](https://floricaconstantine.github.io/TESSERA/reference/normalize.md)
  : Helper function for Moran's I: Center a vector.
- [`plot_neighbor_distances()`](https://floricaconstantine.github.io/TESSERA/reference/plot_neighbor_distances.md)
  : Prepare data for the TESSERA method: Choose a distance threshold for
  adjacency matrices.
- [`poisson_loglike()`](https://floricaconstantine.github.io/TESSERA/reference/poisson_loglike.md)
  : Compute the log likelihood of a set of Poisson variables.
- [`precomp_moran()`](https://floricaconstantine.github.io/TESSERA/reference/precomp_moran.md)
  : Precompute coordinate-dependent spatial constants for Moran's I.
- [`prep_Q_lattice_cpp()`](https://floricaconstantine.github.io/TESSERA/reference/prep_Q_lattice_cpp.md)
  : Pre-align Static Sparse Matrices onto a Shared CSC Template for
  Lattice Q Updates.
- [`prep_data()`](https://floricaconstantine.github.io/TESSERA/reference/prep_data.md)
  : Prepare data for the TESSERA method.
- [`prep_synth_data()`](https://floricaconstantine.github.io/TESSERA/reference/prep_synth_data.md)
  : Generate synthetic multi-sample spatial count data
- [`sample_Poisson_lattice()`](https://floricaconstantine.github.io/TESSERA/reference/sample_Poisson_lattice.md)
  : Simulate spatial counts from a known lattice structure
- [`sample_Poisson_spNNGP()`](https://floricaconstantine.github.io/TESSERA/reference/sample_Poisson_spNNGP.md)
  : Simulate spatial counts from a known spNNGP structure
- [`select_Wald_threshold()`](https://floricaconstantine.github.io/TESSERA/reference/select_Wald_threshold.md)
  : Optimal threshold selection for empirical null estimation
- [`select_Wald_threshold_old()`](https://floricaconstantine.github.io/TESSERA/reference/select_Wald_threshold_old.md)
  : Optimal threshold selection for empirical null estimation
- [`sparse_dist()`](https://floricaconstantine.github.io/TESSERA/reference/sparse_dist.md)
  : Get k nearest neighbors/distances given coordinates.
- [`sparse_dist_LT()`](https://floricaconstantine.github.io/TESSERA/reference/sparse_dist_LT.md)
  : Get k nearest neighbors/distances given coordinates.
- [`sparse_dist_LT_cpp()`](https://floricaconstantine.github.io/TESSERA/reference/sparse_dist_LT_cpp.md)
  : Fast C++ Lower-Triangular Nearest Neighbor Distance Matrix.
- [`sparse_dist_cpp()`](https://floricaconstantine.github.io/TESSERA/reference/sparse_dist_cpp.md)
  : Fast C++ Upper-Triangular Nearest Neighbor Distance Matrix.
- [`summarize_TESSERA()`](https://floricaconstantine.github.io/TESSERA/reference/summarize_TESSERA.md)
  : Summarizes the results of the TESSERA algorithms for a single gene.
- [`tessera_example_data`](https://floricaconstantine.github.io/TESSERA/reference/tessera_example_data.md)
  : Synthetic TESSERA Example Data Files
- [`variables_from_list()`](https://floricaconstantine.github.io/TESSERA/reference/variables_from_list.md)
  : Utility function to instantiate variables from elements in a named
  list. E.g., call this function on list(a=1, b=2) would result in
  variables a and b with values 1 and 2, respectively, in the calling
  frame/environment.
