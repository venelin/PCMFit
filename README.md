
<!-- README.md is generated from README.Rmd. Please edit that file -->
PCMFit
======

The goal of PCMFit is to provide a generic fast maximum likelihood fit and model selection of phylogenetic comparative models (PCMs) to a phylogenetic tree and multivariate trait data at its tips. Supports Gaussian and mixed Gaussian phylogenetic models (MGPM) over all types of trees (including non-ultrametric trees and polytomies). Supports non-existing traits or missing measurements for some of the traits on some of the species. Supports specifying or inferring measurement error for the trait values associated with each tip of the tree. The Gaussian phylogenetic models include various parametrizations of Brownian motion (BM) and Ornstein-Uhlenbeck (OU) multivariate branching processes. The mixed Gaussian models represent models with shifts in the model parameters as well as the type of model at points of the tree. Each shift-point is described as a couple of a shift-node and associated type of model (e.g. OU or BM) driving the trait evolution from the beginning of the branch leading to the shift-node toward the shift-node and its descendants until reaching a tip or another shift-point. The function PCMFit is used to fit a given PCM or a MGPM for a given tree with specified shift-points. The function PCMFitMixed is used to fit an ensemble of possible MGPMs over a tree for which the shift-points are unknown. This function can perform model selection of the best MGPM for a given tree and data according to an information loss estimate such as the Akaike information criterion (AIC).

Installation
------------

### Github

Currently the package can be installed from github using the command:

``` r
devtools::install_github("venelin/PCMFit")
```

### CRAN

Release of the package to CRAN is planned after release of the first stable and documented version.

Documentation
-------------

The writing of the user guide for this R-package is currently in progress. Please, contact the author for assistance, in case you need to use PCMFit.
