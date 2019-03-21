
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Travis build status](https://travis-ci.org/venelin/PCMFit.svg?branch=master)](https://travis-ci.org/venelin/PCMFit) [![Coverage status](https://codecov.io/gh/venelin/PCMFit/branch/master/graph/badge.svg)](https://codecov.io/github/venelin/PCMFit?branch=master) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/PCMFit?color=blue)](https://cran.r-project.org/package=PCMFit) [![Downloads](http://cranlogs.r-pkg.org/badges/PCMFit?color=blue)](https://cran.r-project.org/package=PCMFit)

PCMFit: Statistical inference of phylogenetic comparative models
================================================================

The goal of PCMFit is to provide a generic tool for inference and selection of phylogenetic comparative models (PCMs). Currently, the package implements Gaussian and mixed Gaussian phylogenetic models (MGPM) over all tree types (including non-ultrametric and polytomic trees). The package supports non-existing traits or missing measurements for some of the traits on some of the species. The package supports specifying measurement error associated with each tip of the tree or inferring a measurement error parameter for a group of tips. The Gaussian phylogenetic models include various parametrizations of Brownian motion (BM) and Ornstein-Uhlenbeck (OU) multivariate branching processes. The mixed Gaussian models represent models with shifts in the model parameters as well as the type of model at points of the tree. Each shift-point is described as a pair of a shift-node and associated type of model (e.g. OU or BM) driving the trait evolution from the beginning of the branch leading to the shift-node toward the shift-node and its descendants until reaching a tip or another shift-point. The function PCMFit is used to fit a given PCM or a MGPM for a given tree with specified shift-points. The function PCMFitMixed is used to fit an ensemble of possible MGPMs over a tree for which the shift-points are unknown. This function can perform model selection of the best MGPM for a given tree and data according to an information loss function such as the Akaike information criterion (AIC). The package has been thoroughly tested and applied to real data in the related research article entitled "Automatic Generation of Evolutionary Hypotheses using Mixed Gaussian Phylogenetic Models" (currently in review). Currently, the package is available from <https://github.com/venelin/PCMFit>. The web-page <https://venelin.github.io/PCMFit/> provides access to documentation and related resources.

<!--An early version of this article is available in Chapter 7 of the doctoral thesis available at <https://doi.org/10.3929/ethz-b-000315296>. -->
Installation
------------

### Github

Currently the package can be installed from github using the command:

``` r
devtools::install_github("venelin/PCMFit")
```

### CRAN

Publishing the package on CRAN is planned after release of the first stable and documented version.

Resources
---------

**Note:** *The writing of the vignettes and [help articles](https://venelin.github.io/PCMFit/reference/index.html) for this package is in progress. Please, contact the author for assistance, in case you need to use PCMFit. Thanks for your understanding.*

The user guides and technical reference for the library are available on the [PCMFit web-page](https://venelin.github.io/PCMFit/):

-   The [Getting started](https://venelin.github.io/PCMFit/articles/pcmfit.html) guide introduces mixed Gaussian phylogenetic models (MGPMs) and provides an example how to use the function `PCMFit()` to infer such models on a given tree and trait data.
-   The [Inferring an MGPM with Unknown Shifts](https://venelin.github.io/PCMFit/articles/pcmfitmixed.html) guide shows how to use the function `PCMFitMixed()` to select the best MGPM for a given tree and trait data, based on an information loss function such as the Akaike information criterion (AIC) *(in preparation)*.
-   The [Performing Parametric Bootstrap of an MGPM](https://venelin.github.io/PCMFit/articles/parambootstrap.html) guide shows how to simulate and perform MGPM inference on parametric bootstrap datasets in order to assess the uncertainty of a given MGPM *(in preparation)*.

The research article "Automatic Generation of Evolutionary Hypotheses using Mixed Gaussian Phylogenetic Models" provides a general introduction to MGPMs and reports a real data example and a simulation based comparison of MGPMs versus other implementations of phylogenetic comparative models with shifts. The article is currently undergoing peer review for a publication.

The PCMFit source code is located in the [PCMFit github repository](https://github.com/venelin/PCMFit).

PCMFit builds on top of a stack of three tools enabling fast likelihood calculation and simulation of MGPMs:

-   The R-package [PCMBase](https://venelin.github.io/PCMBase/) implements the specification, likelihood calculation and simulation of MGPMs (Mitov et al. 2018).
-   The auxiliary package [PCMBaseCpp](https://github.com/venelin/PCMBaseCpp) provides a fast C++ implementation of the likelihood calculation as described in (Mitov et al. 2018).
-   PCMBaseCpp relies on the C++ library [SPLITT](https://github.com/venelin/SPLITT) implementing fast traversal of phylogenetic trees (Mitov and Stadler 2018).

Feature requests, bugs, etc can be reported in the [PCMFit issues list](https://github.com/venelin/PCMFit/issues).

Citing PCMFit
=============

To give credit to the PCMFit package in a publication, please cite one of the following:

-   Mitov, V., Bartoszek, K., & Stadler, T. (2018). Automatic Generation of Evolutionary Hypotheses using Mixed Gaussian Phylogenetic Models. Chapter 7, Doctoral Thesis No. [25428](https://doi.org/10.3929/ethz-b-000315296) ETH Zurich.
-   Mitov, V., Bartoszek, K., Asimomitis, G., & Stadler, T. (2018, September 24). Fast likelihood evaluation for multivariate phylogenetic comparative methods: the PCMBase R package. arXiv.org. <https://arxiv.org/abs/1809.09014>.
-   Mitov, V., & Stadler, T. (2018). Parallel Likelihood Calculation for Phylogenetic Comparative Models: the SPLITT C++ Library. bioRxiv, 235739. <http://doi.org/10.1101/235739>

Used R-packages
===============

The PCMFit R-package uses the following 3rd party R-packages:

-   For tree processing in R: ape v5.2 (<span class="citeproc-not-found" data-reference-id="R-ape">**???**</span>), data.table v1.12.0 (<span class="citeproc-not-found" data-reference-id="R-data.table">**???**</span>), PCMBase v1.2.9 (<span class="citeproc-not-found" data-reference-id="R-PCMBase">**???**</span>);
-   For algebraic manipulation: expm v0.999.3 (<span class="citeproc-not-found" data-reference-id="R-expm">**???**</span>), mvtnorm v1.0.10 (<span class="citeproc-not-found" data-reference-id="R-mvtnorm">**???**</span>);
-   For plotting: ggtree v1.14.6 (<span class="citeproc-not-found" data-reference-id="R-ggtree">**???**</span>), ggplot2 v3.1.0 (<span class="citeproc-not-found" data-reference-id="R-ggplot2">**???**</span>);
-   For unit-testing: testthat v2.0.1 (<span class="citeproc-not-found" data-reference-id="R-testthat">**???**</span>);
-   For documentation and web-site generation: roxygen2 v6.1.1 (<span class="citeproc-not-found" data-reference-id="R-roxygen2">**???**</span>), pkgdown v1.3.0 (<span class="citeproc-not-found" data-reference-id="R-pkgdown">**???**</span>);

References
==========

Mitov, Venelin, and Tanja Stadler. 2018. “Parallel likelihood calculation for phylogenetic comparative models: The SPLITT C++ library.” *Methods in Ecology and Evolution*, December, 2041–210X.13136.

Mitov, Venelin, Krzysztof Bartoszek, Georgios Asimomitis, and Tanja Stadler. 2018. “Fast likelihood evaluation for multivariate phylogenetic comparative methods: the PCMBase R package.” *arXiv.org*, September, arXiv:1809.09014. <http://arxiv.org/abs/1809.09014>.
