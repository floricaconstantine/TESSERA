# TESSERA: Tool for Estimating Spatial and Sample-level Effects via Regression Analysis

**TESSERA** is a method based on a spatial generalized linear model for
count data where independent components or disconnected segments
function as individual samples. This methodology allows for the
estimation of a common set of fixed effects across multiple samples
while simultaneously accounting for the unique, independent spatial
correlation structure inherent to each sample. By modeling these spatial
dependencies, **TESSERA** provides a rigorous statistical foundation for
inference compared to models that treat observations as independent or
those that fail to account for the spatial context of the data.

While this method was developed primarily to facilitate multi-sample
gene differential expression analysis in spatial transcriptomics, the
underlying statistical approach is broadly applicable to any count-based
spatial data structured into discrete units. Because the method
estimates a shared set of fixed effects across all samples, it enables
the use of hypothesis testing to address complex comparative questions.
For example, identifying significant differences between distinct
experimental conditions would be statistically infeasible if each sample
were analyzed in isolation.

## Installation

**TESSERA** is currently in the process of submission to Bioconductor.
Until the official release, please install the development version from
[GitHub](https://github.com/floricaconstantine/TESSERA) using the
following code:

``` r
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("floricaconstantine/TESSERA")
```

After acceptance, the release version of **TESSERA** will be available
on [Bioconductor](https://bioconductor.org/packages/TESSERA) and can be
installed with the following code:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("TESSERA")
```

## Documentation

A detailed vignette demonstrating the application of **TESSERA** to
multi-sample spatial transcriptomics data is available as a [rendered
HTML
article](https://floricaconstantine.github.io/TESSERA/articles/vignette.html)).

The raw source code for the vignette can be found in the
[`vignettes/`](https://floricaconstantine.github.io/TESSERA/vignettes/vignette.Rmd)
folder.

## Citation

If you use **TESSERA** in your research, please cite the following
preprint:

Constantine, F., Laszik, Z., Dudoit, S., & Purdom, E. (2026). Unlocking
Multi-Sample Differential Expression for Spatial Transcriptomics Data
with TESSERA. *bioRxiv*.

You can also use the following BibTeX entry:

``` bibtex
@article{constantine2026tessera,
  title={Unlocking Multi-Sample Differential Expression for Spatial Transcriptomics Data with TESSERA},
  author={Constantine, Florica and Laszik, Zoltan and Dudoit, Sandrine and Purdom, Elizabeth},
  journal={bioRxiv},
  year={2026},
}
```
