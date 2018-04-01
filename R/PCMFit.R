# Copyright 2018 Venelin Mitov
#
# This file is part of PCMFit
#
# PCMFit is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PCMFit is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PCMFit.  If not, see <http://www.gnu.org/licenses/>.


#' Numerical lower bound
#' @param model a PCM object
#' @return a PCM object of the same S3 classes as model. Calling
#' \code{\link{PCMBase::PCMSetOrGetVecParams}} on this object returns a lower
#' bound for that can be used, e.g. in a call to \code{\link{optim}}
#' @examples
#' model <- PCM("BM3", 3)
#' PCMLowerBound(model)
#' @export
PCMLowerBound <- function(model, ...) {
  UseMethod("PCMLowerBound", model)
}

#' @export
PCMLowerBound.PCM <- function(model, lowerBoundValue = -10, lowerBoundValuePositiveDiag = 0, ...) {
  if(lowerBoundValuePositiveDiag < 0 ) {
    stop("ERR:04000:PCMFit:PCMFit.R:PCMLowerBound.PCM:: lowerBoundValuePositiveDiag should be non-negative.")
  }
  # a vector with the actual parameters excluding repeated and fixed values
  par <- double(PCMNumParams(model))
  # we set the default upper bound value, but entries corresponding to diagonal
  # elements in Choleski factor upper triangular matrix parameters should be set
  # to lowerBoundValuePositiveDiag
  par[] <- lowerBoundValue

  # all parameters unrolled in a vector including repeated parameter values as
  # well as fixed parameter values
  fullParamVector <- PCMGetVecParamsFull(model)
  maxFullParam <- max(fullParamVector, na.rm = TRUE)
  if(!is.finite(maxFullParam)) {
    maxFullParam <- as.double(1)
  }

  # a tricky way to insert values in the model parameters that are certainly not
  # among the fixed non-countable parameter values. We will use match of these
  # values to assign the either value lowerBoundValuePositiveDiag to the right entries in the
  # model parameter matrices.
  parMask <- maxFullParam + (1:PCMNumParams(model))

  # set the values that match unique positions in par.
  PCMSetOrGetVecParams(model, parMask)
  # Find Choleski factors of positive definite matrices. Such parameters need to
  # have positive diagonal elements, i.e. lowerBoundValuePositiveDiag.
  specParams <- attr(model, "specParams", exact = TRUE)
  for(name in names(specParams)) {
    if(specParams[[name]]$type[1] %in% c("matrix", "gmatrix") &&
       length(specParams[[name]]$type) >= 3 &&
       specParams[[name]]$type[3] == "positive.diag" ) {
      if(specParams[[name]]$type[1] == "gmatrix") {
        mi <- match(diag(model[[name]]), parMask)
        par[unique(mi)] <- lowerBoundValuePositiveDiag
      } else if(specParams[[name]]$type[1] == "matrix") {
        # model[[name]] is a k x k x R array
        R <- PCMNumRegimes(model)
        for(r in 1:R) {
          mi <- match(diag(model[[name]][,,r]), parMask)
          par[unique(mi)] <- lowerBoundValuePositiveDiag
        }
      }
    }
  }
  PCMSetOrGetVecParams(model, par)
  model
}

#' @export
PCMLowerBound.MRG <- function(model, X0 = NULL, Sigmae_x = NULL, ...) {
  model <- NextMethod()
  specParams <- attr(model, "specParams", exact = TRUE)

  if(!is.null(specParams$X0) && !is.null(X0)) {
    model$X0 <- X0
  }
  for(name in names(specParams)) {
    if(specParams[[name]]$type[1]=="model") {
      model[[name]] <- PCMLowerBound(model[[name]], X0 = X0, Sigmae_x = Sigmae_x, ...)
    }
  }
  if(!is.null(specParams$Sigmae_x) && !is.null(Sigmae_x)) {
    model$Sigmae_x <- Sigmae_x
  }
  model
}

#' Numerical upper bound
#' @param model a PC
#' M object
#' @return a PCM object of the same S3 classes as model. Calling
#' \code{\link{PCMBase::PCMSetOrGetVecParams}} on this object returns an upper
#' bound for that can be used, e.g. in a call to \code{\link{optim}}
#' #' model <- PCM("BM3", 3)
#' PCMLowerBound(model)
#' @export
PCMUpperBound <- function(model, ...) {
  UseMethod("PCMUpperBound", model)
}

#' @export
PCMUpperBound.PCM <- function(model, upperBoundValue = 10, upperBoundValuePositiveDiag = 10, ...) {
  if(upperBoundValuePositiveDiag <= 0 ) {
    stop("ERR:04010:PCMFit:PCMFit.R:PCMUpperBound.PCM:: upperBoundValuePositiveDiag should be positive.")
  }
  # a vector with the actual parameters excluding repeated and fixed values
  par <- double(PCMNumParams(model))
  # we set the default upper bound value, but entries corresponding to diagonal
  # elements in Choleski factor upper triangular matrix parameters should be set
  # to upperBoundValuePositiveDiag
  par[] <- upperBoundValue

  # all parameters unrolled in a vector including repeated parameter values as
  # well as fixed parameter values
  fullParamVector <- PCMGetVecParamsFull(model)
  maxFullParam <- max(fullParamVector, na.rm = TRUE)
  if(!is.finite(maxFullParam)) {
    maxFullParam <- as.double(1)
  }

  # a tricky way to insert values in the model parameters that are certainly not
  # among the fixed non-countable parameter values. We will use match of these
  # values to assign the either value upperBoundValuePositiveDiag to the right entries in the
  # model parameter matrices.
  parMask <- maxFullParam + (1:PCMNumParams(model))

  # set the values that match unique positions in par.
  PCMSetOrGetVecParams(model, parMask)
  # Find Choleski factors of positive definite matrices. Such parameters need to
  # have positive diagonal elements, i.e. upperBoundValuePositiveDiag.
  specParams <- attr(model, "specParams", exact = TRUE)
  for(name in names(specParams)) {
    if(specParams[[name]]$type[1] %in% c("matrix", "gmatrix") &&
       length(specParams[[name]]$type) >= 3 &&
       specParams[[name]]$type[3] == "positive.diag" ) {
      if(specParams[[name]]$type[1] == "gmatrix") {
        mi <- match(diag(model[[name]]), parMask)
        par[unique(mi)] <- upperBoundValuePositiveDiag
      } else if(specParams[[name]]$type[1] == "matrix") {
        # model[[name]] is a k x k x R array
        R <- PCMNumRegimes(model)
        for(r in 1:R) {
          mi <- match(diag(model[[name]][,,r]), parMask)
          par[unique(mi)] <- upperBoundValuePositiveDiag
        }
      }
    }
  }
  PCMSetOrGetVecParams(model, par)
  model
}

#' @export
PCMUpperBound.MRG <- function(model, X0 = NULL, Sigmae_x = NULL, ...) {
  model <- NextMethod()
  specParams <- attr(model, "specParams", exact = TRUE)

  if(!is.null(specParams$X0) && !is.null(X0)) {
    model$X0 <- X0
  }
  for(name in names(specParams)) {
    if(specParams[[name]]$type[1]=="model") {
      model[[name]] <- PCMUpperBound(model[[name]], X0 = X0, Sigmae_x = Sigmae_x, ...)
    }
  }
  if(!is.null(specParams$Sigmae_x) && !is.null(Sigmae_x)) {
    model$Sigmae_x <- Sigmae_x
  }
  model
}


#' Fitting a PCM model to a given tree and data
#' @inheritParams PCMBase::PCMLik
#' @return an object of class PCMFit
#' @export
PCMFit <- function(X, tree, model, ...) {
  UseMethod("PCMFit", model)
}

#' @importFrom OptimMCMC runOptimAndMCMC configOptimAndMCMC
#' @importFrom PCMBase PCMCreateLikelihood PCMInfo PCMNumParams
#' @export
PCMFit.PCM <- function(
  X, tree, model, metaI = PCMInfo(X, tree, model),
  lik = NULL, prior = NULL, input.data = NULL, config = NULL,
  argsPCMLowerBound = NULL,
  argsPCMUpperBound = NULL,
  argsPCMSetOrGetVecParams = NULL,
  argsConfigOptimAndMCMC = NULL,
  verbose = FALSE) {

  if(is.null(lik)) {
    lik <- PCMCreateLikelihood(X, tree, model, metaI)
  }
  if(is.null(config)) {
    # lowerVecParams <- double(PCMNumParams(model))
    lowerModel <- do.call(PCMLowerBound, c(list(model = model), argsPCMLowerBound))
    #do.call(PCMSetOrGetVecParams, c(list(model = lowerModel, vecParams = lowerVecParams, set = FALSE), argsPCMSetOrGetVecParams))
    #PCMSetOrGetVecParams(model = lowerModel, vecParams = lowerVecParams, set = FALSE)
    lowerVecParams <- do.call(PCMGetVecParams, c(list(model = lowerModel),
                                                 argsPCMSetOrGetVecParams))

    #upperVecParams <- double(PCMNumParams(model))
    upperModel <- do.call(PCMUpperBound, c(list(model = model), argsPCMUpperBound))
    #do.call(PCMSetOrGetVecParams, c(list(model = upperModel, vecParams = upperVecParams, set = FALSE), argsPCMSetOrGetVecParams))
    #PCMSetOrGetVecParams(model = upperModel, vecParams = upperVecParams, set = FALSE)
    upperVecParams <- do.call(PCMGetVecParams, c(list(model = upperModel), argsPCMSetOrGetVecParams))

    config <- do.call(configOptimAndMCMC, c(list(lik = lik, parLower = lowerVecParams, parUpper = upperVecParams), argsConfigOptimAndMCMC) )
  }

  res <- runOptimAndMCMC(lik, prior, input.data, config, verbose)

  res$X <- X
  res$tree <- tree
  res$modelInit <- model
  res$lik <- lik
  res$config <- config
  res$lowerVecParams <- lowerVecParams
  res$lowerModel <- lowerModel
  res$upperVecParams <- upperVecParams
  res$upperModel <- upperModel

  if(!is.null(res$Optim)) {
    par <- res$Optim$par
    res$modelOptim <- model
    #do.call(PCMSetOrGetVecParams, c(list(model = res$modelOptim, vecParams = par),
     #                               argsPCMSetOrGetVecParams))
    PCMSetOrGetVecParams(model = res$modelOptim, vecParams = par)
    res$logLikOptim <- res$Optim$value
  } else {
    res$modelOptim <- NULL
    res$logLikOptim <- as.double(PCMOptions()$PCMBase.Value.NA)
  }

  class(res) <- c("PCMFit", class(res))
  res
}

#' @export
is.PCMFit <- function(object) {
  inherits(object, "PCMFit")
}

#' @importFrom PCMBase PCMOptions PCMTreeNumTips
#' @export
logLik.PCMFit <- function(object, ...) {
  if(!is.PCMFit(object)) {
    stop("object must inherit from class PCMFit.")
  }

  value <- as.double(PCMOptions()$PCMBase.Value.NA)

  if(!is.null(object$modelOptim)) {
    value <- object$logLikOptim
  }

  attr(value, "df") <- PCMNumParams(object$modelInit, countRegimeChanges = TRUE, countModelTypes = TRUE)
  attr(value, "nobs") <- PCMTreeNumTips(object$tree)

  value
}

#' @importFrom PCMBase PCMGetVecParamsRegimesAndModels
#' @export
coef.PCMFit <- function(object, ...) {
  if(!is.PCMFit(object)) {
    stop("object must inherit from class PCMFit.")
  }

  model <- object$modelInit

  if(!is.null(object$modelOptim)) {
    model <- object$modelOptim
  }

  par <- PCMGetVecParamsRegimesAndModels(model, object$tree, ...)
  c(par, numParam = PCMNumParams(model))
}

#' Generate the next mapping of a set of model-types to a number of regimes
#' @param mapping a vector with elements from \code{models}.
#' (repetitions are allowed) denoting a current model-to-regime mapping.
#' @param models a vector of PCM model-types or integer indices.
#' @return a vector of the same length as mapping with elements from \code{models}.
#' @export
PCMNextMapping <- function(mapping, models) {
  R <- length(mapping)
  mappingInd <- match(mapping, models)
  if(any(is.na(mappingInd))) {
    stop(paste0("ERR:04020:PCMFit:PCMFit.R:PCMNextMapping:: mapping should have
                length ", R, " and contain only elements among models; mapping = (", toString(mapping),")",
                ", models=(", toString(models), "), mappingInd=(", toString(mappingInd), ")."))
  }


  numModels <- length(models)

  carry <- 1
  for(pos in R:1) {
    if(mappingInd[pos] + carry <= numModels) {
      mappingInd[pos] <- mappingInd[pos] + carry
      carry <- 0
    } else {
      mappingInd[pos] <- 1
    }
  }
  res <- models[mappingInd]
  names(res) <- names(mapping)
  res
}

#' Iterator over combinations with repetions of a given set of models
#' @param mapping a vector of elements from models giving the initial combination
#' @param models a vector of unique elements to choose from when building the
#' combinations.
#' @return an iterator object with S3 class c("imapping", "abstractiter", "iter").
#' Calling repeatedly nextElem on this object iterates over all possible combinations
#' with repetitions of the same length as the argument \code{mapping}.
#' @examples
#' it <- PCMIteratorMapping(c("BM3", "BM3"), c("BM3", "OU3", "JOU3"))
#' nextElem(it)
#' nextElem(it)
#' @export
PCMIteratorMapping <- function(mapping, models) {
  state <- new.env()

  state$initial <- mapping
  state$current <- NULL # at initial state before calling nextEl

  nextEl <- function() {
    if(is.null(state$current)) {
      state$current <- state$initial
      state$current
    } else {
      nextMapping <- PCMNextMapping(state$current, models)
      if(isTRUE(all(state$initial==nextMapping))) {
        stop("StopIteration", call. = FALSE)
      } else {
        state$current <- nextMapping
        nextMapping
      }
    }
  }
  obj <- list(nextElem=nextEl)
  class(obj) <- c('PCMIteratorMapping', 'abstractiter', 'iter')
  obj
}


#' Fit all possible model-type mappings to the regimes in a tree.
#'
#' @description This function performs multiple model fits of mixed regime models
#' (MRG) mapping different model-types (e.g. BM and OU) to the regimes in a tree.
#' The produced list of MRG-fits can be used for model selection based on information
#' criterions such as AIC.
#' @inheritParams PCMBase::PCMLik
#' @param model a PCM model used primarily to extract attributes such as the
#' number of traits (k) and the set of model-types to use for the mappings.
#' @param
#' @return an object of class PCMSelect
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom stats rnorm
#' @importFrom PCMBase PCMTreeNumUniqueRegimes PCMTreeGetStartingNodesRegimes PCMTreeTableAncestors PCMTreeSplitAtNode PCMTreeSetDefaultRegime PCMTreeNumUniqueRegimes PCMGetVecParams PCMSetOrGetVecParams PCMNumParams
#' @export
PCMFitModelMappings <- function(
  X, tree, modelTypes = c("BM3", "OU3"),

  metaIFun = PCMInfo, lik = NULL, prior = NULL, input.data = NULL, config = NULL,

  fitToClades = FALSE, cladeRootNodes = PCMTreeGetStartingNodesRegimes(tree),

  cladeFits = NULL,
  numJitterRootCladeFit = 100, sdJitterRootCladeFit = 0.5,
  numJitterAllCladeFits = 100, sdJitterAllCladeFits = 0.5,

  argsMRG = NULL,
  argsPCMLowerBound = NULL,
  argsPCMUpperBound = NULL,
  argsPCMSetOrGetVecParams = NULL,
  argsConfigOptimAndMCMC = NULL,

  printFitVectorsToConsole = FALSE,
  setAttributes = TRUE,

  doParallel = FALSE,
  verbose = FALSE) {

  if(PCMTreeNumUniqueRegimes(tree) <= 0) {
    stop("ERR:04030:PCMFit:PCMFit.R:PCMFitModelMappings:: The tree does not have any regimes. Check tree$edge.regime.")
  }
  if(length(modelTypes) == 0) {
    stop("ERR:04031:PCMFit:PCMFit.R:PCMFitModelMappings:: modelTypes is empty.")
  }

  listPCMOptions <- PCMOptions()

  N <- PCMTreeNumTips(tree)
  startingNodesRegimes <- PCMTreeGetStartingNodesRegimes(tree)
  # assume that all starting nodes start a new regime, no duplicated regimes in
  # different parts of the tree.
  R <- length(startingNodesRegimes)

  if(ncol(X) != N) {
    stop("ERR:04032:PCMFit:PCMFit.R:PCMFitModelMappings:: X should have the same number of columns as the number of tips in tree.")
  }
  if(!is.null(colnames(X)) && !all(colnames(X) == tree$tip.label)) {
    stop("ERR:04033:PCMFit:PCMFit.R:PCMFitModelMappings:: X should either have no column names or these should be
         identical to tree$tip.label.")
  }

  if(fitToClades) {
    if(verbose) {
      cat("Fitting to clades...\n")
    }

    tableAnc <- PCMTreeTableAncestors(tree)

    `%op%` <- if(doParallel) `%dopar%` else `%do%`

    res <- foreach(n = cladeRootNodes, .packages = (.packages()), .combine='c') %op% {
      # set PCMBase options from parent process: necessary if this is executed by
      # a worker process from a cluster.

      options(listPCMOptions)
      treeSplit <- PCMTreeSplitAtNode(tree, n, tableAnc, X)
      clade <- treeSplit$clade
      Xclade <- treeSplit$Xclade

      PCMTreeSetDefaultRegime(clade, 1)

      if(verbose) {
        cat("Fitting model-types (", toString(modelTypes), ") to clade starting at node ", n, " ...\n")
      }
      res <- list()
      res[[as.character(n)]] <-
        PCMFitModelMappings(X = Xclade, tree = clade, modelTypes = modelTypes,
                            metaIFun = metaIFun, lik = lik, prior = prior, input.data = input.data, config = config,
                            fitToClades = FALSE, cladeRootNodes = NULL,
                            cladeFits = cladeFits,
                            numJitterRootCladeFit = numJitterRootCladeFit, sdJitterRootCladeFit = sdJitterRootCladeFit,
                            numJitterAllCladeFits = numJitterAllCladeFits, sdJitterAllCladeFits = sdJitterAllCladeFits,
                            argsMRG = argsMRG,
                            argsPCMLowerBound = argsPCMLowerBound, argsPCMUpperBound = argsPCMUpperBound,
                            argsPCMSetOrGetVecParams = argsPCMSetOrGetVecParams, argsConfigOptimAndMCMC = argsConfigOptimAndMCMC,
                            printFitVectorsToConsole = printFitVectorsToConsole,
                            setAttributes = FALSE,
                            doParallel = doParallel,
                            verbose = verbose)
      res
    }
  } else {
    `%op%` <- if(doParallel) `%dopar%` else `%do%`
    res <- foreach(
      mpp = PCMIteratorMapping(rep(1, PCMTreeNumUniqueRegimes(tree)), 1:length(modelTypes)),
      .packages = (.packages()) ) %op% {
        # set PCMBase options from parent process: necessary if this is executed by
        # a worker process from a cluster.
        options(listPCMOptions)

        # create an mrg model
        model <- do.call(MRG, c(list(k = nrow(X),
                                     modelTypes = modelTypes,
                                     mapping = mpp),
                                argsMRG))

        argsConfigOptimAndMCMC_mpp <- argsConfigOptimAndMCMC

        if(!is.null(cladeFits)) {
          # load listParInitOptim with the parameter vector from the ML fits to all clades
          # starting at startingNodesRegimes
          specParams <- attr(model, "specParams", exact = TRUE)
          subModelsLoaded <- rep(FALSE, length(mpp))

          for(r in 1:length(startingNodesRegimes)) {
            nr <- startingNodesRegimes[r]
            mr <- mpp[r]
            if(is.null(cladeFits[[as.character(nr)]]) ||
               is.null(cladeFits[[as.character(nr)]][[mr]])) {
              stop("ERR:04034:PCMFit:PCMFit.R:PCMFitModelMappings:: missing fitVector for clade starting at node ", nr, "with mapped model-type ", mr, ", while creating initial parameter vector for startingNodesRegimes=c(", toString(startingNodesRegimes),                    ") and mapping c(", toString(mpp), "). Check the cladeFits argument.")
            }
            if(verbose) {
              cat("Loading cladeFit for nr=", nr,"; mr=", mr, "\n")
            }

            fitCladeModel <-
              do.call(PCMLoadMRGFromFitVector,
                      c(list(fitVector = cladeFits[[as.character(nr)]][[mr]],
                             modelTypes = modelTypes,
                             k = nrow(X)),
                        argsMRG))

            for(name in names(fitCladeModel)) {
              if(verbose) {
                cat("Setting model or global parameter with name=", name, "\n")
              }

              if(name %in% names(specParams) &&
                 specParams[[name]]$type[1] %in% c("gscalar", "gvector", "gmatrix")) {
                if( specParams[[name]]$type[2] != "fixed" ) {
                  # set all global non-fixed parameters to the mean of their best fits for the
                  # clades
                  model[[name]] <- model[[name]] + (fitCladeModel[[name]]/R)
                }
              } else if(specParams[[name]]$type[1] == "model" &&
                        class(fitCladeModel[[name]])[1] == modelTypes[mr]) {
                if(!subModelsLoaded[r]) {
                  model[[as.character(r)]] <- fitCladeModel[[name]]
                  subModelsLoaded[r] <- TRUE
                } else {
                  stop("ERR:04035:PCMFit:PCMFit.R:PCMFitModelMappings:: submodel for regime ", r, " was already loaded from a best clade fit.")
                }
              } else {
                stop("ERR:04036:PCMFit:PCMFit.R:PCMFitModelMappings:: Found a member (", name, ") in fitCladeModel starting from node (", nr, ") and with class '", class(fitCladeModel[[name]]), "' which is neither a global parameter nor a model of the needed type (", modelTypes[mr], ").")
              }
            }
          }

          matParamsFromCladeFits <- matrix(PCMGetVecParams(model), 1, PCMNumParams(model), byrow = TRUE)
          matParamsJitterRootCladeFit <- matParamsJitterAllCladeFits <- NULL

          # if there is more than one clade in the tree and numJitterRootCladeFit > 0
          if( !is.null(model[["2"]]) && numJitterRootCladeFit > 0 ) {
            vecParamIndex <- 1:ncol(matParamsFromCladeFits)
            modelIndexParams <- model
            PCMSetOrGetVecParams(modelIndexParams, vecParamIndex)
            vecParamIndexRootClade <- as.integer(PCMGetVecParams(modelIndexParams[["1"]]))
            matParamsJitterRootCladeFit <-
              matrix(matParamsFromCladeFits[1,], 2*numJitterRootCladeFit, ncol(matParamsFromCladeFits), byrow=TRUE)
            for(j in vecParamIndexRootClade) {
              matParamsJitterRootCladeFit[, j] <- rnorm(2 * numJitterRootCladeFit,
                                                        mean = matParamsFromCladeFits[1, j],
                                                        sd = sdJitterRootCladeFit)
            }
          }

          if( numJitterAllCladeFits > 0 ) {
            matParamsJitterAllCladeFits <-
              matrix(matParamsFromCladeFits[1,], 2*numJitterAllCladeFits, ncol(matParamsFromCladeFits), byrow=TRUE)
            for(j in 1:ncol(matParamsFromCladeFits)) {
              matParamsJitterAllCladeFits[, j] <- rnorm(2*numJitterAllCladeFits,
                                                        mean = matParamsFromCladeFits[1, j],
                                                        sd = sdJitterRootCladeFit)
            }
          }

          if(!is.null(matParamsJitterRootCladeFit) || !is.null(matParamsJitterAllCladeFits)) {
            # need to remove the parameters that go out of the lower-upper bound
            lowerModel <- do.call(PCMLowerBound, c(list(model = model), argsPCMLowerBound))
            lowerVecParams <- do.call(PCMGetVecParams, c(list(model = lowerModel), argsPCMSetOrGetVecParams))

            upperModel <- do.call(PCMUpperBound, c(list(model = model), argsPCMUpperBound))
            upperVecParams <- do.call(PCMGetVecParams, c(list(model = upperModel), argsPCMSetOrGetVecParams))

            if( !is.null(matParamsJitterRootCladeFit) ) {
              matParamsJitterRootCladeFit <-
                do.call(rbind,
                        lapply(1:nrow(matParamsJitterRootCladeFit), function(i) {
                          if(isTRUE(all(matParamsJitterRootCladeFit[i,] >= lowerVecParams)) &&
                             isTRUE(all(matParamsJitterRootCladeFit[i,] <= upperVecParams))
                          ) {
                            matParamsJitterRootCladeFit[i,]
                          } else {
                            NULL
                          }
                        }))
              if(nrow(matParamsJitterRootCladeFit) > numJitterRootCladeFit) {
                matParamsJitterRootCladeFit <- matParamsJitterRootCladeFit[1:numJitterRootCladeFit, ]
              }
              matParamsFromCladeFits <- rbind(matParamsFromCladeFits,
                                              matParamsJitterRootCladeFit)
            }

            if( !is.null(matParamsJitterAllCladeFits) ) {
              matParamsJitterAllCladeFits <-
                do.call(rbind,
                        lapply(1:nrow(matParamsJitterAllCladeFits), function(i) {
                          if(isTRUE(all(matParamsJitterAllCladeFits[i,] >= lowerVecParams)) &&
                             isTRUE(all(matParamsJitterAllCladeFits[i,] <= upperVecParams))
                          ) {
                            matParamsJitterAllCladeFits[i,]
                          } else {
                            NULL
                          }
                        }))
              if(nrow(matParamsJitterAllCladeFits) > numJitterAllCladeFits) {
                matParamsJitterAllCladeFits <- matParamsJitterAllCladeFits[1:numJitterAllCladeFits, ]
              }
              matParamsFromCladeFits <- rbind(matParamsFromCladeFits,
                                              matParamsJitterAllCladeFits)
            }
          }


          if(is.null(argsConfigOptimAndMCMC_mpp)) {
            argsConfigOptimAndMCMC_mpp <- list()
          }

          if( !is.null(argsConfigOptimAndMCMC_mpp[["listParInitOptim"]]) ) {
            if(verbose) {
              cat("Prepending the following parameter vectors to argsConfigOptimAndMCMC[['listParInitOptim']] : ")
              print(matParamsFromCladeFits)
            }
            argsConfigOptimAndMCMC_mpp[["listParInitOptim"]] <-
              c(lapply(1:nrow(matParamsFromCladeFits), function(i) matParamsFromCladeFits[i, ]),
                argsConfigOptimAndMCMC_mpp[["listParInitOptim"]])
          } else {
            if(verbose) {
              cat("Setting listParInitOptim to : ")
              print(matParamsFromCladeFits)
            }
            argsConfigOptimAndMCMC_mpp[["listParInitOptim"]] <-
              lapply(1:nrow(matParamsFromCladeFits), function(i) matParamsFromCladeFits[i, ])
          }

          if( !is.null(argsConfigOptimAndMCMC_mpp[["genInitNumEvals"]]) ) {
            argsConfigOptimAndMCMC_mpp[["genInitNumEvals"]] <- argsConfigOptimAndMCMC_mpp[["genInitNumEvals"]] + nrow(matParamsFromCladeFits)
          } else {
            argsConfigOptimAndMCMC_mpp[["genInitNumEvals"]] <- nrow(matParamsFromCladeFits)
          }
        }

        if(verbose) {
          cat("Performing ML fit on mapping: (", toString(mpp), ") of model-types (", modelTypes, ")\n")
          # cat("Initial model:\n", do.call(paste, c(as.list(capture.output(print(model))), sep = "\n")))
          # cat("names(model): (", toString(names(model)), ").\n")
          # cat("attr(model, 'regimes'): (", toString(attr(model, "regimes")), ").\n")
          # cat("PCMRegimes(model, tree): (", toString(PCMRegimes(model, tree)), ").\n")
        }

        fit <- PCMFit(X = X, tree = tree, model = model, metaI = metaIFun,
                      lik = lik, prior = prior, input.data = input.data,
                      config = config,
                      argsPCMLowerBound = argsPCMLowerBound,
                      argsPCMUpperBound = argsPCMUpperBound,
                      argsPCMSetOrGetVecParams = argsPCMSetOrGetVecParams,
                      argsConfigOptimAndMCMC = argsConfigOptimAndMCMC_mpp,
                      verbose = verbose)

        ll <- logLik(fit)
        vec <- c(coef(fit), logLik = ll, df = attr(ll, "df"), nobs = attr(ll, "nobs"), AIC = AIC(fit))

        if(printFitVectorsToConsole) {
          cat(toString(unname(vec)), "\n")
        }
        vec
      }
    res
  }

  class(res) <- "PCMFitModelMappings"
  if(setAttributes) {
    attr(res, "X") <- X
    attr(res, "tree") <- tree
    attr(res, "modelTypes") <- modelTypes
    attr(res, "metaIFun") <- metaIFun
    attr(res, "lik") <- lik
    attr(res, "prior") <- prior
    attr(res, "input.data") <- input.data
    attr(res, "config") <- config
    attr(res, "fitToClades") <- fitToClades
    attr(res, "cladeRootNodes") <- cladeRootNodes
    attr(res, "cladeFits") <- cladeFits
    attr(res, "numJitterRootCladeFit") <- numJitterRootCladeFit
    attr(res, "sdJitterRootCladeFit") <- sdJitterRootCladeFit
    attr(res, "numJitterAllCladeFits") <- numJitterAllCladeFits
    attr(res, "sdJitterAllCladeFits") <- sdJitterAllCladeFits
    attr(res, "argsMRG") <- argsMRG
    attr(res, "argsPCMLowerBound") <- argsPCMLowerBound
    attr(res, "argsPCMUpperBound") <- argsPCMUpperBound
    attr(res, "argsPCMSetOrGetVecParams") <- argsPCMSetOrGetVecParams
    attr(res, "argsConfigOptimAndMCMC") <- argsConfigOptimAndMCMC
    attr(res, "printFitVectorsToConsole") <- printFitVectorsToConsole
    attr(res, "setAttributes") <- setAttributes
    attr(res, "doParallel") <- doParallel
    attr(res, "verbose") <- verbose
  }

  res
}

#' Load a mixed-regime Gaussian model from a fit vector
#' @param fitVector a numeric vector returned by
#' @return an MRG model
#' @export
PCMLoadMRGFromFitVector <- function(fitVector, modelTypes, k, ...) {
  # the last entries in fitVector are in the following order from left to right:
  # numNumericParams, logLik, df, nobs, AIC;
  # the first elements from 1 to numNumericParams are the actual numeric parameters;
  # the entries that follow are between the numeric parameters and numNumericParams.
  # These must be a pair number 2R, where R is the number of regimes in the tree.
  # The first R of these are the starting nodes for the R regimes (starting from the root-node,
  # which is always the starting node for the first regime); The second R are the indices that
  # map modelTypes to the regimes on the tree.
  last <- length(fitVector)
  AIC <- fitVector[last]
  nobs <- fitVector[last - 1]
  df <- fitVector[last - 2]
  ll <- fitVector[last - 3]
  numNumericParams <- fitVector[last - 4]

  mappingModelsToRegimes <- matrix(as.integer(fitVector[(numNumericParams+1):(last-5)]), nrow = 2, byrow = TRUE)
  model <- MRG(k = k, modelTypes = modelTypes, mapping = mappingModelsToRegimes[2,], ...)

  PCMSetOrGetVecParams(model, fitVector)

  attr(model, "startingNodesRegimes") <- mappingModelsToRegimes[1,]
  attr(model, "ll") <- ll
  attr(model, "df") <- df
  attr(model, "nobs") <- nobs
  attr(model, "AIC") <- AIC
  model
}
