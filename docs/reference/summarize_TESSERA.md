# Summarizes the results of the TESSERA algorithms for a single gene.

Summarizes the results of the TESSERA algorithms for a single gene.

## Usage

``` r
summarize_TESSERA(TESSERAData_obj, TESSERAOutput_obj)
```

## Arguments

- TESSERAData_obj:

  Object containing data. Created by the prepData method. Input to the
  TESSERA algorithms.

- TESSERAOutput_obj:

  Output of the TESSERA algorithms.

## Value

A dataframe with MSE, spatial parameters, and Moran's I values. Key
(General): Each row corresponds to one sample. n_cells is the number of
measurements in the sample. gamma_hat/tau2_hat are spatial parameters in
lattice models. kernel_type, nugget_hat, range_hat, and smoothness_hat
are spatial parameters in th spNNGP models. Key for MSE/means: \_sample
indicates that mean/MSE is computed for a single sample. \_total
indicates that mean/MSE is computed across all samples. counts2 refers
to squared counts. Key for Moran's I values: counts: the observed
counts/inputs. predictions: the predicted counts, equal to theta \*
(library size). residuals: counts - predictions. phi: the estimated
spatial random effects. Xbeta: Covariates X multiplied by estimated
effects beta. eta: phi + X beta. theta: Posterior expectation of
exp(eta). librarysize: Often, the total counts per each cell.

## Author

Florica J Constantine, florica AT berkeley.edu
