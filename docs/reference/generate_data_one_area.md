# Simulate spatial count data for a single area

Generates synthetic spatial transcriptomics data based on a Poisson
Generalized Linear Mixed Model (\\GLMM\\). The function simulates
spatial random effects \\\phi\\ using a \\CAR\\, \\SAR\\, or \\Leroux\\
model, generates covariates \\X\\ based on several stochastic processes,
and samples counts \\z\\ from a Poisson distribution.

## Usage

``` r
generate_data_one_area(
  n_points,
  nb_dist,
  model_type,
  beta_true,
  gamma_true = NULL,
  tau2_true = NULL,
  X_type = "rand_bern",
  ar_gamma = 0.75,
  X = NULL,
  library_size = NULL
)
```

## Arguments

- n_points:

  Integer: Number of observations/points in the area.

- nb_dist:

  Numeric: Distance threshold for determining adjacency. Since
  coordinates are sampled in \\(0, 1)^2\\, a value of 0.03 is typically
  appropriate for \\n \approx 1000\\.

- model_type:

  Character: The spatial model for random effects. Options are "CAR",
  "SAR", or "Leroux".

- beta_true:

  Numeric vector: True fixed effect coefficients \\\beta\\.

- gamma_true:

  Numeric: True spatial correlation parameter \\\gamma\\.

- tau2_true:

  Numeric: True spatial scale parameter \\\tau^2\\.

- X_type:

  Character: Method for generating the design matrix \\X\\. Options
  include:

  - **"rand_bern"**: Bernoulli covariates \\X\_{i, j} \in \\0, 1\\\\.

  - **"rand_unif"**: Uniformly distributed \\X\_{i, j} \sim U(0, 1)\\.

  - **"rand_norm"**: Normally distributed \\X\_{i, j} \sim N(0, 1)\\.

  - **"ar1"**: \\X\\ follows a spatial AR(1) process across sorted
    coordinates.

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
  Defaults to a vector of 1s if not provided.

## Value

A list containing the following components:

- **X**: The generated design matrix of covariates.

- **W**: The sparse adjacency matrix.

- **D**: The diagonal degree matrix.

- **eig_list**: Eigenvalues of the standardized adjacency or Laplacian
  matrix, depending on `model_type`.

- **coords**: A matrix of \\(x, y)\\ coordinates.

- **Q**: The unscaled precision matrix defined by \\W\\, \\D\\, and
  \\\gamma\\.

- **phi_true**: The true sampled spatial random effects \\\phi \sim N(0,
  \tau^2 Q^{-1})\\.

- **eta_true**: The linear predictor \\\eta = X\beta + \phi\\.

- **theta_true**: The relative abundance \\\theta = \exp(\eta)\\.

- **z**: The observed count vector sampled from \\Pois(library\\size \*
  \theta)\\.

## Note

This function requires the `Matrix` package for sparse matrix
operations.

## Author

Florica J Constantine, florica AT berkeley.edu

## Examples

``` r
set.seed(2026)
# Generate data for 1000 cells with a Leroux spatial structure
sim_data <- generate_data_one_area(
  n_points = 1000,
  nb_dist = 0.03,
  model_type = "Leroux",
  beta_true = c(1, 0, -1),
  gamma_true = 0.5,
  tau2_true = 1.0,
  X_type = "rand_bern"
)
```
