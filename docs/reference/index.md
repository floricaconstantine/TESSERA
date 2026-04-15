# Package index

## All functions

- [`BRISC_wrapper()`](BRISC_wrapper.md) : Wrapper to fit sparse NN GP
  model via BRISC.
- [`E_step_Vhat()`](E_step_Vhat.md) : Compute the Covariance Matrix of
  the random effects eta. Part of the E-Step in the EM algorithm.
  Follows Clayton and Kaldor (1987) for the Poisson-Gaussian model.
- [`E_step_Vhat_PLU()`](E_step_Vhat_PLU.md) : Compute the INVERSE
  Covariance Matrix of the random effects eta. Part of the E-Step in the
  EM algorithm. Follows Clayton and Kaldor (1987) for the
  Poisson-Gaussian model.
- [`E_step_etahat()`](E_step_etahat.md) : Compute the Expectation of the
  random effects eta. Part of the E-Step in the EM algorithm. Follows
  Clayton and Kaldor (1987) for the Poisson-Gaussian model.
- [`E_step_etahat_PLU()`](E_step_etahat_PLU.md) : Compute the
  Expectation of the random effects eta. Part of the E-Step in the EM
  algorithm. Follows Clayton and Kaldor (1987) for the Poisson-Gaussian
  model.
- [`E_step_predict()`](E_step_predict.md) : Compute the best fit to the
  data. Associated with the E-Step in the EM algorithm, but not needed
  to run. I.e., this is a utility function. Follows Clayton and
  Kaldor (1987) for the Poisson-Gaussian model.
- [`E_step_thetahat()`](E_step_thetahat.md) : Compute the estimated
  poisson parameters theta. Note that eta = log theta. Associated with
  the E-Step in the EM algorithm, but not needed to run. I.e., this is a
  utility function. Follows Clayton and Kaldor (1987) for the
  Poisson-Gaussian model.
- [`E_step_thetahat_PLU()`](E_step_thetahat_PLU.md) : Compute the
  estimated poisson parameters theta. Note that eta = log theta.
  Associated with the E-Step in the EM algorithm, but not needed to run.
  I.e., this is a utility function. Follows Clayton and Kaldor (1987)
  for the Poisson-Gaussian model.
- [`M_step_BRISC()`](M_step_BRISC.md) : Holding beta constant, use BRISC
  to find the spatial parameters of eta - X beta. Part of the M-Step in
  the EM algorithm. ONLY APPLIES TO THE GAUSSIAN PROCESS MODEL.
- [`M_step_beta()`](M_step_beta.md) : Maximize the expected likelihood
  in beta, holding other variables constant. Part of the M-Step in the
  EM algorithm.
- [`M_step_gamma_CAR()`](M_step_gamma_CAR.md) : Maximize the expected
  likelihood in gamma, holding other variables constant. Part of the
  M-Step in the EM algorithm. ONLY APPLIES TO THE CAR MODEL.
- [`M_step_gamma_Leroux()`](M_step_gamma_Leroux.md) : Maximize the
  expected likelihood in gamma, holding other variables constant. Part
  of the M-Step in the EM algorithm. ONLY APPLIES TO THE LEROUX MODEL.
- [`M_step_gamma_SAR()`](M_step_gamma_SAR.md) : Maximize the expected
  likelihood in gamma, holding other variables constant. Part of the
  M-Step in the EM algorithm. ONLY APPLIES TO THE SAR MODEL.
- [`M_step_tau2()`](M_step_tau2.md) : Maximize the expected likelihood
  in tau^2, holding other variables constant. Part of the M-Step in the
  EM algorithm.
- [`M_step_variogram()`](M_step_variogram.md) : Holding beta constant,
  fit a variogram model to eta - X beta. Part of the M-Step in the EM
  algorithm. Instead of a traditional MLE, which has cubic time
  complexity, we fit a variogram to estimate the kernel parameters,
  which is quadratic time. ONLY APPLIES TO THE GAUSSIAN PROCESS MODEL.
- [`Q_matrix_CAR()`](Q_matrix_CAR.md) : Compute the unscaled precision
  matrix in a CAR model.
- [`Q_matrix_Leroux()`](Q_matrix_Leroux.md) : Compute the unscaled
  precision matrix in a Leroux model.
- [`Q_matrix_SAR()`](Q_matrix_SAR.md) : Compute the unscaled precision
  matrix in a SAR model.
- [`TESSERA_lattice()`](TESSERA_lattice.md) : Fit Multi-Sample Poisson
  Spatial GLMM via ECM Algorithm
- [`TESSERA_spNNGP()`](TESSERA_spNNGP.md) : Fit Multi-Sample Poisson
  Spatial GLMM via spNNGP
- [`calc_moran()`](calc_moran.md) : Fast computation of Moran's I.
- [`checkInputsTESSERA()`](checkInputsTESSERA.md) : Check inputs for the
  TESSERA method. Make sure that the input object has everything needed
  to run without error.
- [`checkInputsTESSERAspNNGP()`](checkInputsTESSERAspNNGP.md) : Check
  inputs for the TESSERA_spNNGP method. Make sure that the input object
  has everything needed to run without error.
- [`distanceCalculate()`](distanceCalculate.md) : Helper function for
  Moran's I.
- [`expected_loglike()`](expected_loglike.md) : Compute the expected log
  likelihood.
- [`fit_scaled_noncentral_chi2()`](fit_scaled_noncentral_chi2.md) : Fit
  a scaled non-central chi^2_1 distribution using BOBYQA
- [`generate_data_one_area()`](generate_data_one_area.md) : Simulate
  spatial count data for a single area
- [`generate_data_one_area_spNNGP()`](generate_data_one_area_spNNGP.md)
  : Simulate spatial count data via spNNGP
- [`inversePrecisionMatrixWald()`](inversePrecisionMatrixWald.md) :
  Helper function for Wald stat covariances.
- [`kernel.exp()`](kernel.exp.md) : Exponential Kernel Function.
- [`kernel.gauss()`](kernel.gauss.md) : Gaussian (Squared Exponential)
  Kernel Function.
- [`kernel.matern()`](kernel.matern.md) : Matern Kernel Function.
- [`kernel.sph()`](kernel.sph.md) : Spherical Kernel Function.
- [`moran_I_nb()`](moran_I_nb.md) : Compute Moran's I using an adjacency
  matrix.
- [`neg_hessian_beta()`](neg_hessian_beta.md) : Compute the Negative
  Hessian for beta.
- [`neg_hessian_gamma_CAR()`](neg_hessian_gamma_CAR.md) : Compute the
  Negative Hessian of gamma. ONLY APPLIES TO THE CAR MODEL.
- [`neg_hessian_gamma_Leroux()`](neg_hessian_gamma_Leroux.md) : Compute
  the Negative Hessian of gamma. ONLY APPLIES TO THE Leroux MODEL.
- [`neg_hessian_gamma_SAR()`](neg_hessian_gamma_SAR.md) : Compute the
  Negative Hessian of gamma. ONLY APPLIES TO THE SAR MODEL.
- [`neg_hessian_tau2()`](neg_hessian_tau2.md) : Compute the Negative
  Hessian for tau^2.
- [`nngp_prec_mat()`](nngp_prec_mat.md) : Form a sparse Nearest-Neighbor
  Gaussian Process Precision Matrix. See
  https://mc-stan.org/users/documentation/case-studies/nngp.html.
- [`normalize()`](normalize.md) : Helper function for Moran's I: Center
  a vector.
- [`poisson_loglike()`](poisson_loglike.md) : Compute the log likelihood
  of a set of Poisson variables.
- [`prepData()`](prepData.md) : Prepare data for the TESSERA method.
- [`prepSynthData()`](prepSynthData.md) : Generate synthetic
  multi-sample spatial count data
- [`sample_Poisson_lattice()`](sample_Poisson_lattice.md) : Simulate
  spatial counts from a known lattice structure
- [`sample_Poisson_spNNGP()`](sample_Poisson_spNNGP.md) : Simulate
  spatial counts from a known spNNGP structure
- [`scaledNonCentralChi2PValues()`](scaledNonCentralChi2PValues.md) :
  Given Wald statistics, compute p-values.
- [`selectWaldStatisticThreshold()`](selectWaldStatisticThreshold.md) :
  Optimal threshold selection for empirical null estimation
- [`sparseDist()`](sparseDist.md) : Get k nearest neighbors/distances
  given coordinates. Taken from:
  https://stackoverflow.com/questions/5560218/computing-sparse-pairwise-distance-matrix-in-r
- [`sparseDist_LT()`](sparseDist_LT.md) : Get k nearest
  neighbors/distances given coordinates. Modified from:
  https://stackoverflow.com/questions/5560218/computing-sparse-pairwise-distance-matrix-in-r
- [`summarizeTESSERAPerformance()`](summarizeTESSERAPerformance.md) :
  Summarizes the results of the TESSERA algorithms for a single gene.
- [`variables_from_list()`](variables_from_list.md) : Utility function
  to instantiate variables from elements in a named list. E.g., call
  this function on list(a=1, b=2) would result in variables a and b with
  values 1 and 2, respectively, in the calling frame/environment.
- [`visualizeNeighborDistances()`](visualizeNeighborDistances.md) :
  Prepare data for the TESSERA method: Choose a distance threshold for
  adjacency matrices.
- [`waldStatisticPValuesThreshold()`](waldStatisticPValuesThreshold.md)
  : Given Wald statistics, compute p-values.
- [`waldTestStastics()`](waldTestStastics.md) : Given the output of the
  TESSERA algorithms and a contrast matrix, compute Wald T-statistics.
