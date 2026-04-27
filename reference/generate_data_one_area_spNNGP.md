# Simulate spatial count data via spNNGP

Generates synthetic spatial transcriptomics data using a Poisson
Generalized Linear Mixed Model (\\GLMM\\) with a Sparse Nearest Neighbor
Gaussian Process (\\spNNGP\\) prior for the random effects. This
function approximates a Gaussian Process by using a specified number of
nearest neighbors to construct a sparse precision matrix.

## Usage

``` r
generate_data_one_area_spNNGP(
  n_points,
  nb_dist,
  cov_type,
  cov_params,
  nngp_k,
  beta_true,
  X_type = "rand_bern",
  ar_gamma = 0.75,
  X = NULL,
  library_size = NULL
)
```

## Arguments

- n_points:

  Integer: Number of points/observations in the area.

- nb_dist:

  Numeric: Distance threshold for determining adjacency (if needed for
  secondary lattice structures). Coordinates are sampled in \\(0,
  1)^2\\.

- cov_type:

  Character: The covariance kernel for the Gaussian process. Options are
  "Exp" (Exponential), "Mat" (Matern), "Gau" (Gaussian), and "Sph"
  (Spherical).

- cov_params:

  Numeric vector: Parameters for the covariance kernel (typically
  sigma^2/sill, tau^2/nugget, phi/range, and optionally nu/smoothness).

- nngp_k:

  Integer: The number of nearest neighbors used in the spNNGP
  approximation.

- beta_true:

  Numeric vector: True fixed effect coefficients \\\beta\\.

- X_type:

  Character: Method for generating the design matrix \\X\\. Options
  include:

  - **"rand_bern"**: Bernoulli-style covariates \\X\_{i, j} \in \\0,
    1\\\\.

  - **"rand_unif"**: Uniformly distributed \\X\_{i, j} \sim U(0, 1)\\.

  - **"rand_norm"**: Normally distributed \\X\_{i, j} \sim N(0, 1)\\.

  - **"ar1"**: \\X\\ follows an AR(1) process across sorted coordinates.

  - **"ar1_bern"**: AR(1) process transformed into binary values.

  - **"intercept"**: A constant column of 1s.

- ar_gamma:

  Numeric: Autocorrelation parameter for \\X\\ (used if `X_type` is
  "ar1" or "ar1_bern").

- X:

  Optional matrix: A pre-defined design matrix to use instead of
  generating one.

- library_size:

  Optional numeric vector: Pre-specified library sizes (offsets).
  Defaults to a vector of 1s.

## Value

A list containing the following components:

- **X**: The generated design matrix of covariates.

- **W**: The sparse adjacency matrix (if calculated).

- **D**: The diagonal degree matrix.

- **eig_list**: Eigenvalues for comparison lattice structures (if
  calculated).

- **coords**: A matrix of \\(x, y)\\ coordinates.

- **Q**: The sparse precision matrix \\Q = (I - A)^T D^{-1} (I - A)\\.

- **Dinv**: The diagonal variance matrix \\D^{-1}\\ from the NNGP
  decomposition.

- **A**: The sparse lower triangular matrix representing neighbor
  weights.

- **phi_true**: The sampled spatial random effects \\\phi \sim N(0,
  Q^{-1})\\.

- **eta_true**: The linear predictor \\\eta = X\beta + \phi\\.

- **theta_true**: The relative abundance \\\theta = \exp(\eta)\\.

- **z**: The observed count vector sampled from \\Pois(library\\size \*
  \theta)\\.

## Note

This function requires the `Matrix` package for handling sparse
precision structures. The spNNGP approach is significantly more
memory-efficient than full Gaussian Processes for large \\n\\.

## Author

Florica J Constantine, florica AT berkeley.edu

## Examples

``` r
set.seed(2026)
# Generate data using a Matern kernel with 20 nearest neighbors
sim_data <- generate_data_one_area_spNNGP(
  n_points = 1000,
  nb_dist = 0.03,
  cov_type = "Mat",
  cov_params = c(1.0, 0.1, 10, 2), # sill, nugget, range, smoothness
  nngp_k = 20,
  beta_true = c(1, 0, -1)
)
```
