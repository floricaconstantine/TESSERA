# Simulate spatial counts from a known lattice structure

Generates synthetic counts from a Poisson lattice model (CAR, SAR, or
Leroux) given a pre-defined adjacency structure and fixed parameters.
This is particularly useful for benchmarking against existing spatial
graphs or performing power analyses on specific tissue architectures.

## Usage

``` r
sample_Poisson_lattice(
  model_type,
  X,
  W,
  D,
  library_size,
  tau2_true,
  gamma_true,
  beta_true
)
```

## Arguments

- model_type:

  Character: The spatial model for random effects. Options are "CAR",
  "SAR", or "Leroux".

- X:

  Matrix: Covariates/design matrix.

- W:

  Matrix: Sparse adjacency matrix representing
  measurement-to-measurement connectivity.

- D:

  Matrix: Diagonal degree matrix (row-sums of `W`).

- library_size:

  Numeric vector: Scaling factors (offsets) for each observation.

- tau2_true:

  Numeric: True spatial scale parameter \\\tau^2\\.

- gamma_true:

  Numeric: True spatial correlation parameter \\\gamma\\.

- beta_true:

  Numeric vector: True fixed effect coefficients \\\beta\\.

## Value

A list containing the following components:

- **z**: The observed count vector sampled from \\Pois(library\\size \*
  \theta)\\.

- **phi_true**: The true sampled spatial random effects \\\phi \sim N(0,
  \tau^2 Q^{-1})\\.

- **eta_true**: The linear predictor \\\eta = X\beta + \phi\\.

- **theta_true**: The relative abundance \\\theta = \exp(\eta)\\.

## Note

This function requires the `Matrix` package for sparse matrix operations
and Cholesky decomposition.

## Author

Florica J Constantine, florica AT berkeley.edu

## Examples

``` r
set.seed(2026)
# Use existing area data to re-sample Poisson counts
tau2 <- 1.0
gamma <- 0.5
beta <- c(1, 0, -1)

# Generate base structure
base_struct <- generate_data_one_area(1000, 0.03, "Leroux", beta, gamma, tau2, "rand_bern")

# Re-sample counts from the same structure
sim_counts <- sample_Poisson_lattice(
  model_type = "Leroux",
  X = base_struct$X,
  W = base_struct$W,
  D = base_struct$D,
  library_size = base_struct$library_size,
  tau2_true = tau2,
  gamma_true = gamma,
  beta_true = beta
)
```
