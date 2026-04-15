# Simulate spatial counts from a known spNNGP structure

Generates synthetic counts from a Poisson spNNGP model given a set of
spatial coordinates, covariance parameters, and fixed effects. This
function is ideal for testing the sensitivity of the spNNGP
approximation across different neighborhood sizes \\k\\ or kernel types
without regenerating the underlying spatial layout.

## Usage

``` r
sample_Poisson_spNNGP(
  cov_type,
  X,
  library_size,
  coords,
  cov_params,
  nngp_k,
  beta_true
)
```

## Arguments

- cov_type:

  Character: The covariance kernel for the Gaussian process. Options are
  "Exp" (Exponential), "Mat" (Matern), "Gau" (Gaussian), and "Sph"
  (Spherical).

- X:

  Matrix: Covariates/design matrix where rows correspond to
  observations.

- library_size:

  Numeric vector: Scaling factors (offsets) for each observation.

- coords:

  Matrix: Spatial coordinates \\(x, y)\\ for each observation. Must have
  the same number of rows as `X`.

- cov_params:

  Numeric vector: Covariance parameters ordered as (nugget, sill, range,
  and optionally smoothness).

- nngp_k:

  Integer: The number of nearest neighbors to use for the spNNGP
  precision matrix construction.

- beta_true:

  Numeric vector: True fixed effect coefficients \\\beta\\.

## Value

A list containing the following components:

- **z**: The observed count vector sampled from \\Pois(library\\size \*
  \theta)\\.

- **phi_true**: The true sampled spatial random effects \\\phi \sim N(0,
  Q^{-1})\\, where \\Q\\ is the sparse NNGP precision matrix.

- **eta_true**: The linear predictor \\\eta = X\beta + \phi\\.

- **theta_true**: The relative abundance \\\theta = \exp(\eta)\\.

## Note

This function requires the `Matrix` package for sparse Cholesky
decomposition and solving the precision structure.

## Author

Florica J Constantine, florica AT berkeley.edu

## Examples

``` r
set.seed(2026)
# Simulation parameters
cov_params <- c(0.1, 1.0, 10, 2) # nugget, sill, range, smoothness
beta <- c(1, 0, -1)

# Generate initial structure
base_data <- generate_data_one_area_spNNGP(
  n_points = 1000,
  nb_dist = 0.03,
  cov_type = "Mat",
  cov_params = cov_params,
  nngp_k = 20,
  beta_true = beta
)

# Re-sample counts using the known spNNGP kernel
sim_counts <- sample_Poisson_spNNGP(
  cov_type = "Mat",
  X = base_data$X,
  library_size = base_data$library_size,
  coords = base_data$coords,
  cov_params = cov_params,
  nngp_k = 20,
  beta_true = beta
)
```
