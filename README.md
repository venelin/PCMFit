
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Travis build status](https://travis-ci.org/venelin/PCMFit.svg?branch=master)](https://travis-ci.org/venelin/PCMFit) [![Coverage status](https://codecov.io/gh/venelin/PCMFit/branch/master/graph/badge.svg)](https://codecov.io/github/venelin/PCMFit?branch=master) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/PCMFit?color=blue)](https://cran.r-project.org/package=PCMFit) [![Downloads](http://cranlogs.r-pkg.org/badges/PCMFit?color=blue)](https://cran.r-project.org/package=PCMFit)

<h1>
PCMFit: Statistical inference of phylogenetic comparative models
</h1>
The goal of PCMFit is to provide a generic tool for inference and selection of phylogenetic comparative models (PCMs). Currently, the package implements Gaussian and mixed Gaussian phylogenetic models (MGPM) over all tree types (including non-ultrametric and polytomic trees). The package supports non-existing traits or missing measurements for some of the traits on some of the species. The package supports specifying measurement error associated with each tip of the tree or inferring a measurement error parameter for a group of tips. The Gaussian phylogenetic models include various parametrizations of Brownian motion (BM) and Ornstein-Uhlenbeck (OU) multivariate branching processes. The mixed Gaussian models represent models with shifts in the model parameters as well as the type of model at points of the tree. Each shift-point is described as a pair of a shift-node and associated type of model (e.g. OU or BM) driving the trait evolution from the beginning of the branch leading to the shift-node toward the shift-node and its descendants until reaching a tip or another shift-point. The function PCMFit is used to fit a given PCM or a MGPM for a given tree with specified shift-points. The function PCMFitMixed is used to fit an ensemble of possible MGPMs over a tree for which the shift-points are unknown. This function can perform model selection of the best MGPM for a given tree and data according to an information loss function such as the Akaike information criterion (AIC). The package has been thoroughly tested and applied to real data in the related research article entitled "Automatic Generation of Evolutionary Hypotheses using Mixed Gaussian Phylogenetic Models" (currently in review). Currently, the package is available from <https://github.com/venelin/PCMFit>. The web-page <https://venelin.github.io/PCMFit/> provides access to documentation and related resources.

<!--An early version of this article is available in Chapter 7 of the doctoral thesis available at <https://doi.org/10.3929/ethz-b-000315296>. -->
Prerequisites
=============

Before installing PCMFit, it is necessary to ensure that several R-packages are installed or can be installed from CRAN. These are listed below:

-   [PCMBase](https://venelin.github.io/PCMBase). The PCMBase package is available on CRAN and should be installed automatically during the installation of PCMFit. If this does not happen, try the command:

``` r
install.packages("PCMBase")
```

-   \[[PCMBaseCpp](https://github.com/venelin/PCMBaseCpp)\]. This package contains Rcpp modules for likelihood calculation of the model types implemented in PCMBase. PCMBaseCpp can be used as a companion and not a substitute of PCMBase. The sole purpose of PCMBaseCpp is to accelerate the likelihood calculation by implementing the most computationally intensive algorithm in C++, which has shown a dramatic speed-up (in the order of 100 times) (Mitov et al. 2018). Hence, installing this package is optional but highly recommended, in particular, if the goal is to infer models with shifts and/or to infer models on trees bigger than 100 tips. Installing PCMBaseCpp requires a C++ compiler to be installed on the system. This installation has been tested on two systems:

    -   On Mac OS X, it has been tested using the default (clang) and the Intel (icpc) compiler.
    -   On Linux (Euler ETH cluster), it has been tested with the following modules loaded (commands in the file `.bashrc`):

``` sh
module load interconnect/ethernet
module load new gcc/6.3.0 open_mpi/1.6.5 r/3.4.0 intel/2017.5
```

Upon validating the availability of a C++ compiler, PCMBaseCpp can be installed using the commands:

``` r
# These two packages are available on CRAN
install.packages("Rcpp")
install.packages("RcppArmadillo")
# At the time of writing, PCMBaseCpp is available only from github. 
devtools::install_github("venelin/PCMBaseCpp")
```

-   other third party dependencies include the packages [`data.table`](https://CRAN.R-project.org/package=data.table), [`foreach`](https://CRAN.R-project.org/package=foreach), [`iterators`](https://CRAN.R-project.org/package=iterators), [`ape`](https://CRAN.R-project.org/package=ape) and [`digest`](https://CRAN.R-project.org/package=digest). These packages should be installed automatically from CRAN when installing PCMFit. If this does not happen, consult the packages' web-pages (links above).

-   An optional but highly recommended dependency. The function `PCMTreePlot` in the package PCMBase uses the R-package ggtree, which is not on CRAN. It is highly recommended to install this package in order to be able to visualize trees with colored parts corresponding to defferent evolutionary regimes. If `ggtree` is not installed, the figures in the vignettes and coding examples cannot be generated. At the time of writing this documentation, `ggtree` can be installed from bioconductor through the following code (if that does not work, check the [ggtree home page](https://guangchuangyu.github.io/software/ggtree/)):

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree", version = "3.8")
```

Installing PCMFit
=================

From Github
-----------

Currently the package can be installed from github using the command:

``` r
devtools::install_github("venelin/PCMFit")
```

From CRAN
---------

Publishing PCMFit on CRAN is planned after release of the first stable and documented version.

Parallel execution
==================

Currently PCMFit implements parallel execution for the inference of mixed Gaussian phylogenetic models with unknown shifts. This is optional but highly recommended. To enable parallel execution, it is necessary to run PCMFit on a computer equipped with a multiple core processor or on multiple node computing cluster. In its current implementation, PCMFit uses the function `%dopar%` from the R-package [`foreach`](https://CRAN.R-project.org/package=foreach) to parallelize the execution of (nested) `foreach` loops. This parallelization has been tested using two parallel backends for the `%dopar%` function:

-   Using the R packages [`doMPI`](https://CRAN.R-project.org/package=doMPI) and [`Rmpi`](https://CRAN.R-project.org/package=Rmpi) on a multiple node cluster with `open_mpi/1.6.5` installed. In particular, MGPM inference has been run using up to 250 cores on the [ETH scientific computing cluster Euler](https://scicomp.ethz.ch/wiki/Euler).
-   Using the R package [`doParallel`](https://CRAN.R-project.org/package=doParallel) on a MacBook Pro (Retina, 15-inch, Late 2013), 2.3 GHz Intel Core i7 processor (4 physical cores, 8 logical cores), running macOS Sierra 10.12.6.

To install the above packages, follow the most recent instructions in their documentation (links to the packages web-pages provided above). Once you have installed the parallel backend of choice, you can paste/edit the following code snippet in the beginning of the R-script for running PCMFit model inference:

``` r
library(PCMBase)
library(PCMBaseCpp)
library(PCMFit)
# other needed packages, e.g. ape, data.table etc...

# extract dataset identifier and possibly other parameters from the command line:
args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
  data_id <- as.integer(args[1])
} else {
  data_id <- 1L
}

# A character string used in filenames for a model inference on a given data:
prefixFiles = paste0("MGPM_A_F_BC2_RR_DATAID_", data_id)


# creating the cluster for this PCMFit run:
if(!exists("cluster") || is.null(cluster)) {
  if(require(doMPI)) {
    # using MPI cluster as distributed node cluster (possibly running on a 
    # cluster of multiple nodes)
    # Get the number of cores. Assume this is run in a batch job.
    p = strtoi(Sys.getenv('LSB_DJOB_NUMPROC'))
    cluster <- startMPIcluster(count = p-1, verbose = TRUE)
    doMPI::registerDoMPI(cluster)
  } else {
    # possibly running on personal computer without mpi installation
    cluster <- parallel::makeCluster(
      parallel::detectCores(logical = TRUE),
      outfile = paste0("log_", prefixFiles, ".txt"))
    doParallel::registerDoParallel(cluster)
  }
}
```

Finally, to tell PCMFit that it should run the inference in parallel, specify the argument `doParallel=TRUE` in calls to the function `PCMFitMixed`. A full example for this is provided in the user guide [Inferring an MGPM with Unknown Shifts](https://venelin.github.io/PCMFit/articles/pcmfitmixed.html).

Resources
=========

A brief historical background and theoretical overview of PCMs can be found in Chapter 1 of (Mitov 2018a). A more general introduction can be found in (Harmon 2018). The research article "Automatic Generation of Evolutionary Hypotheses using Mixed Gaussian Phylogenetic Models" provides a general introduction to MGPMs and reports a real data example and a simulation based comparison of MGPMs versus other implementations of phylogenetic comparative models with shifts. The article is currently undergoing peer review for a publication.

The user guides and technical reference for the library are available on the [PCMFit web-page](https://venelin.github.io/PCMFit/):

**Note:** *The writing of the user gudes and [help articles](https://venelin.github.io/PCMFit/reference/index.html) for this package is in progress. Please, contact the author for assistance, in case you need to use PCMFit right away and need help with the coding examples. Thanks for your understanding.*

-   The [Getting started](https://venelin.github.io/PCMFit/articles/pcmfit.html) guide introduces mixed Gaussian phylogenetic models (MGPMs) and provides an example how to use the function `PCMFit()` to infer such models on a given tree and trait data.
-   The [Inferring an MGPM with Unknown Shifts](https://venelin.github.io/PCMFit/articles/pcmfitmixed.html) guide shows how to use the function `PCMFitMixed()` to select the best MGPM for a given tree and trait data, based on an information loss function such as the Akaike information criterion (AIC) *(in preparation)*.
-   The [Performing Parametric Bootstrap of an MGPM](https://venelin.github.io/PCMFit/articles/parambootstrap.html) guide shows how to simulate and perform MGPM inference on parametric bootstrap datasets in order to assess the uncertainty of a given MGPM *(in preparation)*.

The PCMFit source code is located in the [PCMFit github repository](https://github.com/venelin/PCMFit).

Feature requests, bugs, etc can be reported in the [PCMFit issues list](https://github.com/venelin/PCMFit/issues).

Related tools
=============

PCMFit builds on top of a stack of three tools enabling fast likelihood calculation and simulation of MGPMs:

-   The R-package [PCMBase](https://venelin.github.io/PCMBase/) implements the specification, likelihood calculation and simulation of MGPMs (Mitov et al. 2018).
-   The auxiliary package [PCMBaseCpp](https://github.com/venelin/PCMBaseCpp) provides a fast C++ implementation of the likelihood calculation as described in (Mitov et al. 2018).
-   PCMBaseCpp relies on the C++ library [SPLITT](https://github.com/venelin/SPLITT) implementing fast traversal of phylogenetic trees (Mitov and Stadler 2018).

Citing PCMFit
=============

To acknowledge the PCMFit package in a publication or other presentation, please cite one of the following:

-   Mitov, V. (2018). Phylogenetic Comparative Methods in the Era of Big Data, Doctoral Thesis No. [25428](https://doi.org/10.3929/ethz-b-000315296) ETH Zurich.

        @phdthesis{Mitov:2018hh,
        author = {Mitov, Venelin},
        title = {{Phylogenetic Comparative Methods in the Era of Big Data}},
        school = {ETH Zurich},
        year = {2018},
        address = {Zurich},
        url = {https://doi.org/10.3929/ethz-b-000315296},
        }

    Mitov, V., Bartoszek, K., & Stadler, T. (2018). Automatic Generation of Evolutionary Hypotheses using Mixed Gaussian Phylogenetic Models. Article in Review, preprint available in Chapter 7, Doctoral Thesis No. [25428](https://doi.org/10.3929/ethz-b-000315296) ETH Zurich.

        @article{Mitov:2019agh,
        author = {Mitov, Venelin and Bartoszek, Krzysztof and Stadler, Tanja},
        title = {{Automatic Generation of Evolutionary Hypotheses using Mixed Gaussian Phylogenetic Models}},
        journal = {article in review, preprint in Ch. 7 of Thesis No. 25428, ETH Zurich},
        year = {2019},
        pages = {},
        month = mai,
        url = {https://doi.org/10.3929/ethz-b-000315296}
        }

-   Mitov, V., Bartoszek, K., Asimomitis, G., & Stadler, T. (2018, September 24). Fast likelihood evaluation for multivariate phylogenetic comparative methods: the PCMBase R package. arXiv.org. <https://arxiv.org/abs/1809.09014>.

        @article{Mitov:2018fl,
        author = {Mitov, Venelin and Bartoszek, Krzysztof and Asimomitis, Georgios and Stadler, Tanja},
        title = {{Fast likelihood evaluation for multivariate phylogenetic comparative methods: the PCMBase R package}},
        journal = {arXiv.org},
        year = {2018},
        eprint = {1809.09014},
        eprinttype = {arxiv},
        eprintclass = {q-bio.PE},
        pages = {arXiv:1809.09014},
        month = sep,
        annote = {34 pages, 6 figures}
        }

-   Mitov, V., & Stadler, T. (2018). Parallel Likelihood Calculation for Phylogenetic Comparative Models: the SPLITT C++ Library. Methods in Ecology and Evolution, 2041--210X.13136. <http://doi.org/10.1101/235739>

        @article{Mitov:2018dqa,
        author = {Mitov, Venelin and Stadler, Tanja},
        title = {{Parallel likelihood calculation for phylogenetic comparative models: The SPLITT C++ library}},
        journal = {Methods in Ecology and Evolution},
        year = {2018},
        pages = {2041--210X.13136},
        month = dec
        }

Used software packages
======================

Although, I have been consistent in my effort to update the following list with any new package I have used in developing and testing PCMFit, chances are that I have omitted some of these tools. I apologise to their authors.

The PCMFit R-package uses the following 3rd party R-packages:

-   For tree processing in R: ape v5.3 (Paradis et al. 2019), data.table v1.12.0 (Dowle and Srinivasan 2019), PCMBase v1.2.10 (Mitov 2019a);
-   For specification and manipulation of models in R: PCMBase v1.2.10 (Mitov 2019a), PCMBaseCpp v0.1.4 (Mitov 2019b), SPLITT v1.2.1 (Mitov 2018b);
-   For data processing in R: data.table v1.12.0 (Dowle and Srinivasan 2019);
-   For parallel execution: iterators v1.0.10 (Analytics and Weston 2018), foreach v1.4.4 (Revolution Analytics and Weston, n.d.), doParallel v1.0.14 (Corporation and Weston 2018);
-   For algebraic computation: expm v0.999.4 (Goulet et al. 2019), mvtnorm v1.0.10 (Genz et al. 2019);
-   For plotting: ggtree v1.14.6 (Yu and Lam 2019), ggplot2 v3.1.0 (Wickham et al. 2018), cowplot v0.9.4 (Wilke 2019);
-   For unit-testing: testthat v2.0.1 (Wickham 2018);
-   For documentation and web-site generation: roxygen2 v6.1.1 (Wickham, Danenberg, and Eugster 2018), pkgdown v1.3.0 (Wickham and Hesselberth 2018);

References
==========

Analytics, Revolution, and Steve Weston. 2018. *Iterators: Provides Iterator Construct for R*. <https://CRAN.R-project.org/package=iterators>.

Corporation, Microsoft, and Steve Weston. 2018. *DoParallel: Foreach Parallel Adaptor for the ’Parallel’ Package*. <https://CRAN.R-project.org/package=doParallel>.

Dowle, Matt, and Arun Srinivasan. 2019. *Data.table: Extension of ‘Data.frame‘*. <https://CRAN.R-project.org/package=data.table>.

Genz, Alan, Frank Bretz, Tetsuhisa Miwa, Xuefei Mi, and Torsten Hothorn. 2019. *Mvtnorm: Multivariate Normal and T Distributions*. <https://CRAN.R-project.org/package=mvtnorm>.

Goulet, Vincent, Christophe Dutang, Martin Maechler, David Firth, Marina Shapira, and Michael Stadelmann. 2019. *Expm: Matrix Exponential, Log, ’Etc’*. <https://CRAN.R-project.org/package=expm>.

Harmon, Luke J. 2018. *Phylogenetic Comparative Methods*. Learning from Trees. <https://lukejharmon.github.io/pcm/>.

Mitov, Venelin. 2018a. “Phylogenetic Comparative Methods in the Era of Big Data.” PhD thesis, Zurich: ETH Zurich. <https://doi.org/10.3929/ethz-b-000315296>.

———. 2018b. *SPLITT: A Generic Library for Serial and Parallel Lineage Traversal of Trees*.

———. 2019a. *PCMBase: Simulation and Likelihood Calculation of Phylogenetic Comparative Models*. <https://CRAN.R-project.org/package=PCMBase>.

———. 2019b. *PCMBaseCpp: A C++ Implementation of Parallel Likelihood Calculation for Phylogenetic Comparative Models*.

Mitov, Venelin, and Tanja Stadler. 2018. “Parallel likelihood calculation for phylogenetic comparative models: The SPLITT C++ library.” *Methods in Ecology and Evolution*, December, 2041–210X.13136.

Mitov, Venelin, Krzysztof Bartoszek, Georgios Asimomitis, and Tanja Stadler. 2018. “Fast likelihood evaluation for multivariate phylogenetic comparative methods: the PCMBase R package.” *arXiv.org*, September, arXiv:1809.09014. <http://arxiv.org/abs/1809.09014>.

Paradis, Emmanuel, Simon Blomberg, Ben Bolker, Joseph Brown, Julien Claude, Hoa Sien Cuong, Richard Desper, et al. 2019. *Ape: Analyses of Phylogenetics and Evolution*. <https://CRAN.R-project.org/package=ape>.

Revolution Analytics, and Steve Weston. n.d. *Foreach: Provides Foreach Looping Construct for R*.

Wickham, Hadley. 2018. *Testthat: Unit Testing for R*. <https://CRAN.R-project.org/package=testthat>.

Wickham, Hadley, and Jay Hesselberth. 2018. *Pkgdown: Make Static Html Documentation for a Package*. <https://CRAN.R-project.org/package=pkgdown>.

Wickham, Hadley, Winston Chang, Lionel Henry, Thomas Lin Pedersen, Kohske Takahashi, Claus Wilke, and Kara Woo. 2018. *Ggplot2: Create Elegant Data Visualisations Using the Grammar of Graphics*. <https://CRAN.R-project.org/package=ggplot2>.

Wickham, Hadley, Peter Danenberg, and Manuel Eugster. 2018. *Roxygen2: In-Line Documentation for R*. <https://CRAN.R-project.org/package=roxygen2>.

Wilke, Claus O. 2019. *Cowplot: Streamlined Plot Theme and Plot Annotations for ’Ggplot2’*. <https://CRAN.R-project.org/package=cowplot>.

Yu, Guangchuang, and Tommy Tsan-Yuk Lam. 2019. *Ggtree: An R Package for Visualization and Annotation of Phylogenetic Trees with Their Covariates and Other Associated Data*. <https://guangchuangyu.github.io/software/ggtree>.
