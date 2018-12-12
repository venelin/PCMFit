#' Generate the next mapping of a set of model-types to a number of regimes
#' @param mapping a vector with elements from \code{models}.
#' (repetitions are allowed) denoting a current model-to-regime mapping.
#' @param modelTypes a vector of PCM model-types or integer indices.
#' @return a vector of the same length as mapping with elements from \code{modelTypes}.
#' @export
PCMNextMapping <- function(mapping, modelTypes) {
  R <- length(mapping)
  mappingInd <- match(mapping, modelTypes)
  if(any(is.na(mappingInd))) {
    stop(paste0("ERR:04100:PCMFit:PCMFit.R:PCMNextMapping:: mapping should have
                length ", R,
                " and contain only elements among modelTypes; mapping = (",
                toString(mapping),")",
                ", modelTypes=(", toString(modelTypes), "), mappingInd=(",
                toString(mappingInd), ")."))
  }


  numModels <- length(modelTypes)

  carry <- 1
  for(pos in R:1) {
    if(mappingInd[pos] + carry <= numModels) {
      mappingInd[pos] <- mappingInd[pos] + carry
      carry <- 0
    } else {
      mappingInd[pos] <- 1
    }
  }
  res <- modelTypes[mappingInd]
  names(res) <- names(mapping)
  res
}

#' Iterator over combinations with repetions of a given set of modelTypes
#' @param mapping a vector of elements from modelTypes giving the initial combination
#' @param modelTypes a vector of unique elements to choose from when building the
#' combinations.
#' @return an iterator object with S3 class c("imapping", "abstractiter", "iter").
#' Calling repeatedly nextElem on this object iterates over all possible combinations
#' with repetitions of the same length as the argument \code{mapping}.
#' @examples
#' it <- PCMIteratorMapping(c("BM", "BM"), c("BM", "OU", "JOU"))
#' iterators::nextElem(it)
#' iterators::nextElem(it)
#' @export
PCMIteratorMapping <- function(mapping, modelTypes) {
  state <- new.env()

  state$initial <- mapping
  state$current <- NULL # at initial state before calling nextEl

  nextEl <- function() {
    if(is.null(state$current)) {
      state$current <- state$initial
      state$current
    } else {
      nextMapping <- PCMNextMapping(state$current, modelTypes)
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


#' Generate the next mapping of model-types chosen from regime-specific sets of possible model-types.
#' @param mapping a vector with elements from \code{modelTypes}.
#' (repetitions are allowed) denoting a current model-to-regime mapping.
#' @param modelTypes a vector of unique model-types, e.g. \code{c("BM", "OU")}.
#' @param allowedModelTypesIndices a list of the same length as \code{mapping} with integer vector elements
#' or NULLs. When an element of this list is an integer vector its elements denote unique positions in
#' modelTypes, i.e. the allowed model-types for the regime at that position in mapping.
#' @return a vector of the same length as mapping with elements from \code{modelTypes}.
#' @export
PCMNextMapping2 <- function(mapping, modelTypes, allowedModelTypesIndices) {
  R <- length(mapping)
  numModels <- length(modelTypes)

  if(length(allowedModelTypesIndices) != length(mapping)) {
    stop(paste0("ERR:04130:PCMFit:PCMFit.R:PCMNextMapping2:: mapping and allowedModelTypesIndices should
                be the same length, ", R, "."))
  }

  mappingInd <- match(mapping, modelTypes)
  #cat("mappingInd: ", toString(mappingInd), "\n")
  mappingInd2 <- sapply(1:R, function(pos) {
    if(!is.null(allowedModelTypesIndices[[pos]])) {
      #cat("pos = ", pos, "; allowedModelTypesIndicies[[pos]] = ", toString(allowedModelTypesIndices[[pos]]), "\n")
      #cat("pos = ", pos, "; mappingInd[pos] = ", mappingInd[pos], "\n")
      match(mappingInd[pos], allowedModelTypesIndices[[pos]])
    } else {
      # allowedModelTypesIndices[[pos]] == NULL means that all model types are allowed at this position
      match(mappingInd[pos], 1:numModels)
    }
  })


  if(any(is.na(mappingInd2))) {
    stop(paste0("ERR:04100:PCMFit:PCMFit.R:PCMNextMapping2:: mapping should have
                length ", R, " allowedModelTypesIndices; mapping = (", toString(mapping),")",
                ", modelTypes=(", toString(modelTypes), "), mappingInd2=(", toString(mappingInd2), ")."))
  }


  numsAllowedModelTypes <- sapply(1:R, function(pos) {
    if(is.null(allowedModelTypesIndices[[pos]])) {
      numModels
    } else {
      length(allowedModelTypesIndices[[pos]])
    }
  })


  carry <- 1
  for(pos in R:1) {
    if(mappingInd2[pos] + carry <= numsAllowedModelTypes[pos]) {
      mappingInd2[pos] <- mappingInd2[pos] + carry
      carry <- 0
    } else {
      mappingInd2[pos] <- 1
    }
  }


  #cat("mappingInd2=", toString(mappingInd2), "\n")
  res <- sapply(1:R, function(pos) {
    if(is.null(allowedModelTypesIndices[[pos]])) {
      modelTypes[mappingInd2[pos]]
    } else {
      modelTypes[allowedModelTypesIndices[[pos]][mappingInd2[pos]]]
    }
  })

  #cat("indices=", toString(indices), "\n")
  #res <- modelTypes[indices]
  names(res) <- names(mapping)
  res
}

#' Iterator over combinations with repetions of a given set of modelTypes
#' @param mapping a vector of elements from modelTypes giving the initial combination
#' @param modelTypes a vector of unique elements to choose from when building the
#' combinations.
#' @param allowedModelTypesIndices a list of the same length as \code{mapping} with integer vector elements
#' or NAs. When an element of this list is an integer vector its elements denote unique positions in
#' modelTypes, i.e. the allowed model-types for the regime at that position in mapping.
#' @return an iterator object with S3 class c("imapping", "abstractiter", "iter").
#' Calling repeatedly nextElem on this object iterates over all possible combinations
#' with repetitions of the same length as the argument \code{mapping}.
#' @examples
#' it <- PCMIteratorMapping(c("BM", "BM"), c("BM", "OU", "JOU"))
#' iterators::nextElem(it)
#' iterators::nextElem(it)
#' @export
PCMIteratorMapping2 <- function(mapping, modelTypes, allowedModelTypesIndices) {
  state <- new.env()

  state$initial <- mapping
  state$current <- NULL # at initial state before calling nextEl

  nextEl <- function() {
    if(is.null(state$current)) {
      state$current <- state$initial
      state$current
    } else {
      nextMapping <- PCMNextMapping2(state$current, modelTypes, allowedModelTypesIndices)
      if(isTRUE(all(state$initial==nextMapping))) {
        stop("StopIteration", call. = FALSE)
      } else {
        state$current <- nextMapping
        nextMapping
      }
    }
  }
  obj <- list(nextElem=nextEl)
  class(obj) <- c('PCMIteratorMapping2', 'abstractiter', 'iter')
  obj
}


#' Load a mixed-regime Gaussian model from a fit vector
#' @export
PCMLoadMixedGaussianFromFitVector <- function(fitVector, modelTypes, k,
                                              remapModelTypeIndicesInFitVector = seq_along(modelTypes), ...) {
  # the last entries in fitVector are in the following order from left to right:
  # numNumericParams, logLik, df, nobs, AIC;
  # the first elements from 1 to numNumericParams are the actual numeric parameters;
  # the entries that follow are between the numeric parameters and numNumericParams.
  # These must be a pair number 2R, where R is the number of regimes in the tree.
  # The first R of these are the starting nodes for the R regimes (starting from the root-node,
  # which is always the starting node for the first regime); The second R are the indices that
  # map modelTypes to the regimes on the tree.
  last <- length(fitVector)
  v_aic <- fitVector[last]
  nobs <- fitVector[last - 1]
  df <- fitVector[last - 2]
  ll <- fitVector[last - 3]
  numNumericParams <- fitVector[last - 4]

  mappingModelsToRegimes <- matrix(as.integer(fitVector[(numNumericParams+1):(last-5)]), nrow = 2, byrow = TRUE)

  model <- MixedGaussian(k = k, modelTypes = modelTypes,
                         mapping = remapModelTypeIndicesInFitVector[mappingModelsToRegimes[2, ]], ...)

  PCMParamLoadOrStore(model, fitVector, offset = 0, load = TRUE)

  attr(model, "startingNodesRegimes") <- mappingModelsToRegimes[1,]
  attr(model, "ll") <- ll
  attr(model, "df") <- df
  attr(model, "nobs") <- nobs
  attr(model, "aic") <- v_aic
  model
}

#' Match a model mapping vector against a vector of model types
#' @param modelMapping a character or integer vector
#' @param modelTypes a character vector with valid model-class-names.
#' @return a character vector with elements from modelTypes or stops with an error
#' @export
MatchModelMapping <- function(modelMapping, modelTypes) {
  if(is.character(modelMapping)) {
    m <- match(modelMapping, modelTypes)
    if( any(is.na(m)) ) {
      stop("ERR:04111:PCMFit:PCMFit.R:MatchModelMapping:: some of the models in modelMapping could not be matched against model-types in tableFits (", toString(modelMapping[which(is.na(m))]), ")")
    }
  } else if(is.integer(modelMapping)) {
    m <- match(modelMapping, 1:length(modelTypes))
    if( any(is.na(m)) ) {
      stop("ERR:04112:PCMFit:PCMFit.R:MatchModelMapping:: some of the integer models in modelMapping could not be matched against model-indices in tableFits (", toString(modelMapping[which(is.na(m))]), ")")
    }
  } else {
    stop("ERR:04113:PCMFit:PCMFit.R:MatchModelMapping:: modelMapping should be character or integer (",
         toString(modelMapping), ")")
  }
  unname(modelTypes[m])
}

#' @importFrom PCMBase PCMTreeGetStartingNodesRegimes PCMTreeGetLabels
#' @importFrom digest digest
#' @export
HashCodes <- function(tree, modelTypes, startingNodesRegimesLabels, modelMapping) {
  orderSNRL <- order(as.integer(startingNodesRegimesLabels))
  list(
    hashCodeTree = digest(PCMTreeToString(tree), serialize = FALSE),
    hashCodeStartingNodesRegimesLabels = digest(
      toString(startingNodesRegimesLabels[orderSNRL]), serialize = FALSE),
    hashCodeMapping = digest(
      toString(MatchModelMapping(modelMapping[orderSNRL], modelTypes)), serialize = FALSE)
  )
}

#' Lookup a fit vector for a given tree and model mapping in a data.table of
#' previously run fits.
#'
#' @param tree a phylo object
#' @param modelTypes character vector
#' @param modelMapping an integer or character vector to by matched against modelTypes
#' @param tableFits a data.table having at least the following columns:
#' \itemize{
#' \item{hashCodeTree}{an MD5 key column of type character-vector}
#' \item{hashCodeStartingNodesRegimesLabels}{an MD5 key column of type character-vector representing the hash-code of
#' \code{PCMTreeGetLabels(tree)[PCMTreeGetStartingNodesRegimes(tree)]}.}
#' \item{hashCodeMapping}{an MD5 key column of type character-vector}}
#' @return the corresponding fit-vector to the given tree and model mapping or
#' if no such entry is found, issues an error.
#' @importFrom digest digest
#' @export
LookupFit <- function(tree, modelTypes, modelMapping, tableFits,
                      hashCodes = HashCodes(tree = tree,
                                            modelTypes = modelTypes,
                                            startingNodesRegimesLabels =
                                              PCMTreeGetLabels(tree)[PCMTreeGetStartingNodesRegimes(tree)],
                                            modelMapping = modelMapping )) {
  tableFits[hashCodes, , mult="first", nomatch=0]
}

#' Compose a MixedGaussian model
#' @importFrom PCMBase PCMTreeExtractClade is.Fixed is.PCM PCMParamLoadOrStore PCMParamGetShortVector
#' @export
ComposeMixedGaussianFromFits <- function(
  tree, startingNodesRegimes, modelTypes, k, R, mapping,
  argsMixedGaussian,
  tableFits, modelTypesInTableFits,
  tableAncestors = NULL, verbose = FALSE) {

  if(verbose) {
    cat("Composing a MixedGaussian model:\n")
    cat("startingNodesRegimes=c(", toString(startingNodesRegimes), ")\n")
    cat("mapping=c(", toString(mapping), ")\n")
  }

  # tableFits can be supplied from a previous fit with different modelTypes which partially
  # overlap with modelTypes.
  # In this case, we need to remap the model indices in fitVector to the new modelTypes.
  remapModelTypeIndicesInFitVector <- match(modelTypesInTableFits, modelTypes)
  if(verbose) {
    cat("remapModelTypeIndicesInFitVector=c(", toString(remapModelTypeIndicesInFitVector), ");")
    cat("NAs indicate that some of modelTypesInTableFits are not present in modelTypes.\n")
  }

  # create a MixedGaussian model
  model <- do.call(MixedGaussian, c(list(k = k, modelTypes = modelTypes, mapping = mapping), argsMixedGaussian))

  # load listParInitOptim with the parameter vector from the ML fits to all clades
  # starting at startingNodesRegimes
  spec <- attr(model, "spec", exact = TRUE)
  subModelsLoaded <- rep(FALSE, length(mapping))

  for(r in 1:length(startingNodesRegimes)) {
    nr <- startingNodesRegimes[r]
    mr <- modelTypes[mapping[r]]
    tree_nr <- PCMTreeExtractClade(tree, nr, tableAncestors = tableAncestors)
    PCMTreeSetDefaultRegime(tree_nr, 1)
    fit <- LookupFit(tree = tree_nr, modelTypes = modelTypes, modelMapping = mr, tableFits = tableFits)

    if(nrow(fit) == 0) {
      stop("ERR:04120:PCMFit:PCMFit.R:ComposeMixedGaussianFromCladeFits:: no entry in tableFits for the given tree and modelMapping.")
    }

    fitVector <- fit$fitVector[[1]]

    if(verbose) {
      cat("Loading cladeFit for nr=", nr,"; mr=", mr, "\n")
    }

    fitModel <-
      do.call(PCMLoadMixedGaussianFromFitVector,
              c(list(fitVector = fitVector,
                     modelTypes = modelTypes,
                     k = k,
                     remapModelTypeIndicesInFitVector = remapModelTypeIndicesInFitVector),
                argsMixedGaussian))

    for(name in names(fitModel)) {
      if(verbose) {
        cat("Setting model or global parameter with name=", name, "\n")
      }

      if( name %in% names(model) &&
          is.Global(model[[name]]) &&
          is.Global(fitModel[[name]]) &&
          !is.Fixed(model[[name]]) ) {

        # set all global non-fixed parameters to the mean of their best fits for the
        # clades
        vecCurrent <- PCMParamGetShortVector(model[[name]], k = k, R = 1)
        vecFitModel <- PCMParamGetShortVector(fitModel[[name]], k = k, R = 1)

        if(! length(vecCurrent) == length(vecFitModel) ) {
          stop(paste0("ERR:04121:PCMFit:PCMFit.R:ComposeMixedGaussianFromCladeFits:: a global parameter ", name,
                      " in a fit-model has a different short-vector length from the to-be fit model length(vecToFit)=",
                      length(vecCurrent), ", length(vecFitModel)=", length(vecFitModel), "; class model to fit = ", mr,
                      "; class fitted model = ", toString(class(fitModel[[name]]))))
        }
        vecCurrent <- vecCurrent + vecFitModel/R
        PCMParamLoadOrStore(model[[name]], vecCurrent, 0, k, 1, load = TRUE)
        #model[[name]][] <- model[[name]][] + (fitModel[[name]][]/R)

      } else if(is.PCM(spec[[name]]) && class(fitModel[[name]])[1] == mr) {
        if(!subModelsLoaded[r]) {
          model[[as.character(r)]] <- fitModel[[name]]
          subModelsLoaded[r] <- TRUE
        } else {
          stop("ERR:04122:PCMFit:PCMFit.R:ComposeMixedGaussianFromCladeFits:: submodel for regime ", r, " was already loaded from a best clade fit.")
        }
      } else {
        stop("ERR:04123:PCMFit:PCMFit.R:ComposeMixedGaussianFromCladeFits:: Found a member (", name, ") in fitModel starting from node (", nr, ") and with class '", class(fitModel[[name]]), "' which is neither a global parameter nor a model of the needed type (", mr, ").")
      }
    }
  }
  model
}

#' @importFrom PCMBase PCMParamGetShortVector PCMParamLoadOrStore
#' @importFrom stats rnorm
#' @export
AdaptArgsConfigOptimAndMCMC <- function(
  model,

  argsPCMParamLowerLimit,
  argsPCMParamUpperLimit,

  argsPCMParamLoadOrStore,
  argsConfigOptimAndMCMC,
  numJitterRootRegimeFit, sdJitterRootRegimeFit,
  numJitterAllRegimeFits, sdJitterAllRegimeFits,
  verbose = FALSE
) {

  matParamsFromTableFits <- matrix(PCMParamGetShortVector(model), 1, PCMParamCount(model), byrow = TRUE)
  matParamsJitterRootCladeFit <- matParamsJitterAllCladeFits <- NULL

  # if there is more than one clade in the tree and numJitterRootRegimeFit > 0
  if( !is.null(model[["2"]]) && numJitterRootRegimeFit > 0 ) {
    vecParamIndex <- 1:ncol(matParamsFromTableFits)
    modelIndexParams <- model
    PCMParamLoadOrStore(modelIndexParams, vecParamIndex, offset = 0, load = TRUE)
    vecParamIndexRootClade <- as.integer(PCMParamGetShortVector(modelIndexParams[["1"]]))
    matParamsJitterRootCladeFit <-
      matrix(matParamsFromTableFits[1,], 2*numJitterRootRegimeFit, ncol(matParamsFromTableFits), byrow=TRUE)
    for(j in vecParamIndexRootClade) {
      matParamsJitterRootCladeFit[, j] <- rnorm(2 * numJitterRootRegimeFit,
                                                mean = matParamsFromTableFits[1, j],
                                                sd = sdJitterRootRegimeFit)
    }
  }

  if( numJitterAllRegimeFits > 0 ) {
    matParamsJitterAllCladeFits <-
      matrix(matParamsFromTableFits[1,], 2*numJitterAllRegimeFits, ncol(matParamsFromTableFits), byrow=TRUE)
    for(j in 1:ncol(matParamsFromTableFits)) {
      matParamsJitterAllCladeFits[, j] <- rnorm(2*numJitterAllRegimeFits,
                                                mean = matParamsFromTableFits[1, j],
                                                sd = sdJitterAllRegimeFits)
    }
  }

  if(!is.null(matParamsJitterRootCladeFit) || !is.null(matParamsJitterAllCladeFits)) {
    # need to remove the parameters that go out of the lower-upper bound
    lowerModel <- do.call(PCMParamLowerLimit, c(list(model), argsPCMParamLowerLimit))
    lowerVecParams <- PCMParamGetShortVector(lowerModel)

    upperModel <- do.call(PCMParamUpperLimit, c(list(model), argsPCMParamUpperLimit))
    upperVecParams <- PCMParamGetShortVector(upperModel)

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
      if(!is.null(matParamsJitterRootCladeFit) &&
         nrow(matParamsJitterRootCladeFit) > numJitterRootRegimeFit) {
        matParamsJitterRootCladeFit <- matParamsJitterRootCladeFit[1:numJitterRootRegimeFit, ]
      }
      matParamsFromTableFits <- rbind(matParamsFromTableFits,
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
      if(!is.null(matParamsJitterAllCladeFits) &&
         nrow(matParamsJitterAllCladeFits) > numJitterAllRegimeFits) {
        matParamsJitterAllCladeFits <- matParamsJitterAllCladeFits[1:numJitterAllRegimeFits, ]
      }
      matParamsFromTableFits <- rbind(matParamsFromTableFits,
                                      matParamsJitterAllCladeFits)
    }
  }


  if(is.null(argsConfigOptimAndMCMC)) {
    argsConfigOptimAndMCMC <- list()
  }

  if( !is.null(argsConfigOptimAndMCMC[["listParInitOptim"]]) ) {
    if(verbose) {
      cat("Prepending the ", nrow(matParamsFromTableFits), " parameter vectors to argsConfigOptimAndMCMC[['listParInitOptim']] \n ")
      #print(matParamsFromTableFits)
    }
    argsConfigOptimAndMCMC[["listParInitOptim"]] <-
      c(lapply(1:nrow(matParamsFromTableFits), function(i) matParamsFromTableFits[i, ]),
        argsConfigOptimAndMCMC[["listParInitOptim"]])
  } else {
    if(verbose) {
      cat("Setting ", nrow(matParamsFromTableFits), " param vectors into listParInitOptim \n ")
      #print(matParamsFromTableFits)
    }
    argsConfigOptimAndMCMC[["listParInitOptim"]] <-
      lapply(1:nrow(matParamsFromTableFits), function(i) matParamsFromTableFits[i, ])
  }

  if( !is.null(argsConfigOptimAndMCMC[["genInitNumEvals"]]) ) {
    argsConfigOptimAndMCMC[["genInitNumEvals"]] <- argsConfigOptimAndMCMC[["genInitNumEvals"]] + nrow(matParamsFromTableFits)
  } else {
    argsConfigOptimAndMCMC[["genInitNumEvals"]] <- nrow(matParamsFromTableFits)
  }

  argsConfigOptimAndMCMC
}


#' Fit regime-assignments to (sub-)trees in a tree with different assigned model types to each regime.
#'
#' @description This function performs multiple model fits of mixed regime models
#' (MixedGaussian) mapping different model-types (e.g. BM and OU) to different regimes in a tree and
#' testing different regime assignments to the branches in the tree (TODO describe algorithm).
#' @importFrom foreach foreach when %do% %dopar% %:%
#' @importFrom data.table data.table rbindlist is.data.table setkey :=
#' @importFrom PCMBase PCMTreeSetLabels PCMTreeSetDefaultRegime PCMTreeEvalNestedEDxOnTree PCMTreeNumTips PCMTreeListCladePartitions PCMTreeToString MixedGaussian PCMOptions PCMTreeTableAncestors PCMTreeSplitAtNode PCMTreeSetRegimes PCMGetVecParamsRegimesAndModels
#' @importFrom stats logLik coef AIC
#' @return an S3 object of class PCMFitModelMappings.
#'
#' @export
PCMFitModelMappings <- function(
  X, tree, modelTypes,
  SE = matrix(0.0, nrow(X), PCMTreeNumTips(tree)),

  generatePCMModelsFun = NULL,
  metaIFun = PCMInfo, positiveValueGuard = Inf,

  lik = NULL, prior = NULL, input.data = NULL, config = NULL,

  fitMappingsPrev = NULL,
  tableFitsPrev = fitMappingsPrev$tableFits,
  modelTypesInTableFitsPrev = NULL,

  skipFitWhenFoundInTableFits = TRUE,

  prefixFiles = "fits_",

  maxCladePartitionLevel = 8, maxNumNodesPerCladePartition = 1, minCladeSizes = 25,

  listPCMOptions = PCMOptions(),

  argsMixedGaussian = NULL,
  argsPCMParamLowerLimit = NULL,
  argsPCMParamUpperLimit = NULL,
  argsPCMParamLoadOrStore = NULL,

  argsConfigOptimAndMCMC1 = NULL,
  argsConfigOptimAndMCMC2 = NULL,

  numJitterRootRegimeFit = 100, sdJitterRootRegimeFit = 0.5,
  numJitterAllRegimeFits = 100, sdJitterAllRegimeFits = 0.5,

  printFitVectorsToConsole = FALSE,
  setAttributes = TRUE,

  doParallel = FALSE,

  verbose = TRUE,
  verbosePCMFit = FALSE,
  verboseComposeMixedGaussianFromFits = FALSE,
  verboseAdaptArgsConfigOptimAndMCMC = FALSE
) {

  treeOriginal <- tree

  treeEDExpression = "tree"
  PCMTreeSetLabels(tree)
  PCMTreeSetDefaultRegime(tree, 1)

  colnames(X) <- as.character(1:PCMTreeNumTips(tree))

  tableFits <- tableFitsPrev
  if(is.null(modelTypesInTableFitsPrev)) {
    # assume that tableFitsPrev has the same modelTypes unless the user specifies otherwise
    modelTypesInTableFits <- modelTypes
  }


  if(!is.null(fitMappingsPrev)) {
    tableFits <- tableFitsPrev <- fitMappingsPrev$tableFits
    if(!identical(modelTypes, fitMappingsPrev$arguments$modelTypes)) {

      # this should remap the model-type indices in the fit vectors, show
      # table fits is correctly converted.
      tableFits <- RetrieveFittedModelsFromFitVectors(
        fitMappings = fitMappingsPrev, tableFits = tableFitsPrev,
        modelTypesNew = modelTypes)

      modelTypesInTableFits <- modelTypes
    }
  }
  if(!is.data.table(tableFits)) {
    if(verbose) {
      cat("Initiating tableFits...\n")
    }
    tableFits = data.table(treeEDExpression = character(),
                           hashCodeTree = character(),
                           hashCodeStartingNodesRegimesLabels = character(),
                           hashCodeMapping = character(),
                           startingNodesRegimesLabels = list(),
                           mapping = list(),
                           fitVector = list(),
                           logLik = double(),
                           df = integer(),
                           nobs = integer(),
                           aic = double(),
                           duplicated = logical())
    modelTypesInTableFits <- modelTypes
  }
  setkey(tableFits, hashCodeTree,hashCodeStartingNodesRegimesLabels,hashCodeMapping )

  if(length(minCladeSizes) < maxCladePartitionLevel) {
    minCladeSizes <- c(rep(as.integer(NA), maxCladePartitionLevel - length(minCladeSizes)),
                       minCladeSizes)
  }
  MIN_CLADE_SIZE <- min(minCladeSizes, na.rm = TRUE)

  mainLoopHistory <- list()

  if(PCMTreeNumTips(tree) > MIN_CLADE_SIZE) {

    preorderTree <- PCMTreePreorder(tree)
    tableAncestors <- PCMTreeTableAncestors(tree, preorder = preorderTree)


    # 1. (fitsToClades) Perform a fit of each model-type to each clade

    # start from a list containing the trivial partition into one clade equal to
    # the whole tree
    cladeRoots <- c(PCMTreeNumTips(tree) + 1,
                    unlist(PCMTreeListCladePartitions(
                      tree = tree,
                      nNodes = 1,
                      minCladeSize = MIN_CLADE_SIZE,
                      tableAncestors = tableAncestors)))
    if(verbose) {
      cat("Step 1: perform a fit on", length(cladeRoots), " clades x", length(modelTypes), " model-types...\n")
    }

    envCombineTaskResults1 <- new.env()
    envCombineTaskResults1$ncalls <- 0

    SaveCurrentResults <- function(..., filePrefix) {
      status <- try(
        save(..., file = paste0("CurrentResults_", filePrefix, ".RData")),
       silent = TRUE)
      if(class(status) == 'try-error') {
        warning(paste0("An error occurred while saving tableFits to file ",
                       paste0("CurrentResults_", filePrefix, ".RData"), " :", status))
      }
    }

    # combine fits from parallel tasks into a data.table and saves this data.table to
    # an .RData file
    CombineTaskResults <- function(..., filePrefix, envNCalls) {
      envNCalls$ncalls <- envNCalls$ncalls + 1
      data <- rbindlist(
        lapply(
          list(...),
          function(dt) {
            if(is.data.table(dt)) {
              dt
            } else {
              cat("ERROR in foreach task: no data.table returned: \n", toString(dt), "\n")
              NULL
            }
          }))
      save(data, file = paste0(filePrefix, envNCalls$ncalls, ".RData"))
      #data <- data[duplicated == FALSE]
      data
    }

    CleanTemporaryFitFiles <- function(filePrefix) {
      status <- try(
        {
          files <- list.files(".", pattern=paste0("^", filePrefix, ".*.RData"))
          # cat("Cleaning temporary data-files:")
          # print(files)
          if(length(files) > 0) {
            do.call(file.remove, as.list(files))
          }

        }, silent = TRUE)
      if(class(status) == 'try-error') {
        warning(paste0("An error occurred while cleaning up temporary data-files: ", status))
      }
    }


    nodeLabelsTree <- PCMTreeGetLabels(tree)
    `%op%` <- if(doParallel) `%dopar%` else `%do%`

    fitsToClades <- foreach(
      clRoot = cladeRoots,
      clEDExpression = paste0("E(", treeEDExpression, ",", nodeLabelsTree[cladeRoots], ")"),
      .combine = function(...) rbindlist(list(...)),
      .multicombine=TRUE,
      .inorder = FALSE,
      .packages = (.packages()) ) %:%

      foreach(
        modelMapping = 1:length(modelTypes),
        .combine = function(...) {
          CombineTaskResults(
            ...,
            filePrefix = paste0(prefixFiles, "_clades_"),
            envNCalls = envCombineTaskResults1)
        },
        .multicombine = TRUE,
        .inorder = FALSE,
        .errorhandling = "pass",
        .packages = (.packages()) ) %op% {


          # recreate tableAncestors under a different name, to avoid exporting a potentially big
          # tableAncestors
          if(!exists("tableAncestorsTree", .GlobalEnv)) {
            if(verbose) {
              cat("Creating tableAncestorsTree in .GlobalEnv during clade-fits\n")
            }
            assign("tableAncestorsTree", PCMTreeTableAncestors(tree, preorderTree), .GlobalEnv)
          }

          if(!exists("generatedPCMModels", .GlobalEnv) && !is.null(generatePCMModelsFun)) {
            if(verbose) {
              cat("Calling generatePCMModelsFun()...\n")
            }
            # this should generate PCMParentClasses and PCMSpecify functions
            # for all models in modelTypes
            generatePCMModelsFun()
            assign("generatedPCMModels", TRUE, .GlobalEnv)
          }

          # don't want the names in modelMapping
          modelMapping <- unname(modelMapping)

          options(listPCMOptions)
          treeSplit <- PCMTreeSplitAtNode(tree, clRoot, tableAncestorsTree, X)
          clade <- treeSplit$clade
          Xclade <- treeSplit$Xclade

          PCMTreeSetDefaultRegime(clade, 1)

          hashCodes <- HashCodes(
            tree = clade,
            modelTypes = modelTypes,
            startingNodesRegimesLabels = PCMTreeGetLabels(clade)[PCMTreeNumTips(clade) + 1],
            modelMapping = modelMapping)

          fit <- LookupFit(tableFits = tableFits, hashCodes = hashCodes)

          if(nrow(fit) == 1 && skipFitWhenFoundInTableFits) {
            dt.row <- fit
            if(verbose) {
              cat("  Found a fit in tableFits for ", modelTypes[modelMapping], " on clade starting at node ", clRoot, " ...\n")
            }
            if(printFitVectorsToConsole) {
              cat(dt.row$treeEDExpression[[1]], ":", dt.row$hashCodeTree[[1]], ":", dt.row$hashCodeStartingNodesRegimes[[1]], ":", dt.row$hashCodeMapping[[1]], ":", toString(unname(dt.row$fitVector[[1]])), "\n", sep="")
            }
            dt.row[, duplicated:=TRUE]

          } else {
            if(verbose) {
              cat("  Fitting ", modelTypes[modelMapping], " to clade starting at node ", clRoot, " ...\n")
            }

            # create an MixedGaussian model
            model <- do.call(MixedGaussian, c(list(k = nrow(Xclade),
                                         modelTypes = modelTypes,
                                         mapping = modelMapping),
                                    argsMixedGaussian))

            fit <- PCMFit(X = Xclade, tree = clade, model = model, SE = SE,
                          metaI = metaIFun, positiveValueGuard = positiveValueGuard,

                          lik = lik, prior = prior, input.data = input.data,
                          config = config,
                          argsPCMParamLowerLimit = argsPCMParamLowerLimit,
                          argsPCMParamUpperLimit = argsPCMParamUpperLimit,
                          argsPCMParamLoadOrStore = argsPCMParamLoadOrStore,
                          argsConfigOptimAndMCMC = argsConfigOptimAndMCMC1,
                          verbose = verbosePCMFit)

            ll <- unname(logLik(fit))
            v_aic = unname(AIC(fit))
            vec <- c(coef(fit), logLik = ll, df = attr(ll, "df"), nobs = attr(ll, "nobs"), aic = v_aic)
            dt.row <- data.table(
              treeEDExpression = clEDExpression,
              hashCodeTree = hashCodes$hashCodeTree,
              hashCodeStartingNodesRegimesLabels = hashCodes$hashCodeStartingNodesRegimesLabels,
              hashCodeMapping = hashCodes$hashCodeMapping,
              # character names in global tree
              startingNodesRegimesLabels = list(PCMTreeGetLabels(clade)[PCMTreeNumTips(clade) + 1]),
              mapping = list(MatchModelMapping(modelMapping, modelTypes)),
              fitVector = list(unname(vec)),
              logLik = ll,
              df = attr(ll, "df"),
              nobs = attr(ll, "nobs"),
              aic = v_aic,
              duplicated = FALSE)

          }
          if(printFitVectorsToConsole) {
            cat(dt.row$treeEDExpression[[1]], ":", dt.row$hashCodeTree[[1]], ":",
                dt.row$hashCodeStartingNodesRegimes[[1]], ":", dt.row$hashCodeMapping[[1]], ":",
                toString(unname(dt.row$fitVector[[1]])), "\n", sep="")
          }
          dt.row
        }

    setkey(fitsToClades, hashCodeTree,hashCodeStartingNodesRegimesLabels,hashCodeMapping)

    # update tableFits with the entries in fitsToClades which were not already there.
    if(nrow(tableFits)>0) {
      tableFits <- rbindlist(list(tableFits, fitsToClades))
    } else {
      tableFits <- fitsToClades
    }
    tableFits <- tableFits[duplicated == FALSE]
    setkey(tableFits, hashCodeTree,hashCodeStartingNodesRegimesLabels,hashCodeMapping)

    SaveCurrentResults(tableFits, filePrefix = prefixFiles)
    CleanTemporaryFitFiles(filePrefix = prefixFiles)


    # 2. Perform fits to clade-partitions with different model mappings

    # we need these variables throughout this step
    treeRootInt <- PCMTreeNumTips(tree) + 1
    treeRootLabel <- nodeLabelsTree[treeRootInt]
    hashCodeEntireTree <- tableFits[
      sapply(startingNodesRegimesLabels, length) == 1 &
        sapply(startingNodesRegimesLabels, function(s) match(treeRootInt, s, nomatch = 0L)) == 1,
      hashCodeTree[1]]

    # 2.1 identify the best clade-partition and mapping currently in tableFits. This should be the
    # best fit of one of the model types to the whole tree

    queuePartitionRoots <- fitsToClades[
      list(hashCodeEntireTree), {
        iMinAIC <- which.min(aic)
        list(
          level = 1L,
          node = treeRootInt,
          partitionParentNode = treeRootInt,
          hashCodeBestPartitionInitial = hashCodeStartingNodesRegimesLabels[iMinAIC],
          hashCodeBestMappingInitial = hashCodeMapping[iMinAIC],
          hashCodeBestPartitionLevel = hashCodeStartingNodesRegimesLabels[iMinAIC],
          hashCodeBestMappingLevel = hashCodeMapping[iMinAIC]
        )
      }]
    #setkey(queuePartitionRoots, level, node)
    # index of the head row in queuePartitionRoots
    headQPR <- 1L

    envCombineTaskResults2 <- new.env()
    envCombineTaskResults2$ncalls <- 0

    # Main loop:
    while(headQPR <= nrow(queuePartitionRoots)) {
      # 2.2. pop the first partition root from the queue
      partitionRootNode <- queuePartitionRoots[headQPR, node]
      partitionRootLabel <- nodeLabelsTree[partitionRootNode]
      partitionRootLevel <- queuePartitionRoots[headQPR, level]

      hashCodeBestPartition <- queuePartitionRoots[headQPR, hashCodeBestPartitionLevel]

      bestPartition <- tableFits[queuePartitionRoots[headQPR,
                                                     list(hashCodeEntireTree,
                                                          hashCodeBestPartitionLevel,
                                                          hashCodeBestMappingLevel)],
                                 startingNodesRegimesLabels[[1]]]
      bestMapping <- tableFits[queuePartitionRoots[headQPR,
                                                   list(hashCodeEntireTree,
                                                        hashCodeBestPartitionLevel,
                                                        hashCodeBestMappingLevel)],
                               mapping[[1]]]

      bestAIC <- tableFits[queuePartitionRoots[headQPR,
                                               list(hashCodeEntireTree,
                                                    hashCodeBestPartitionLevel,
                                                    hashCodeBestMappingLevel)],
                           aic[[1]]]

      # prepare an entry in the mainLoopHistory for this iteration
      mainLoopHistoryEntry <- list(
        headQPR = headQPR,
        lengthQPR = nrow(queuePartitionRoots),
        headQPR_PartitionRootNode = partitionRootNode,
        headQPR_PartitionRootLabel = partitionRootLabel,
        headQPR_PartitionRootLevel = partitionRootLevel,
        headQPR_Partition = bestPartition,
        headQPR_MappingIdx = match(bestMapping, modelTypes),
        headQPR_Mapping = bestMapping,
        headQPR_AIC = bestAIC)

      if(verbose) {
        cat("Step 2.2: headQPR=", headQPR, ": bestPartition/mapping/AIC: ", toString(bestPartition), " / ", toString(match(bestMapping, modelTypes)), "(", toString(bestMapping), ") / ", bestAIC, " ; \n")
        cat("  Remaining queue:\n")

        print(cbind(
          queuePartitionRoots[-(1:headQPR), list(level, node, partitionParentNode)],
          tableFits[queuePartitionRoots[-(1:headQPR),
                                        list(hashCodeEntireTree,
                                             hashCodeBestPartitionInitial,
                                             hashCodeBestMappingInitial)],
                    list(partInit = startingNodesRegimesLabels)],
          tableFits[queuePartitionRoots[-(1:headQPR),
                                        list(hashCodeEntireTree,
                                             hashCodeBestPartitionLevel,
                                             hashCodeBestMappingLevel)],
                    list(partCurrent = startingNodesRegimesLabels)]))

        cat("Performing clade-partitioning at partitionRootLabel=", partitionRootLabel, "; partitionRootLevel=", partitionRootLevel, "\n")
      }

      # advance the head to the next entry in the queue
      headQPR <- headQPR + 1

      # 2.3. Extract the tree
      edExpression <- paste0("E(tree,", partitionRootLabel, ")")
      for(label in bestPartition) {
        if(tableAncestors[label, partitionRootLabel] > 0) {
          edExpression <- paste0("D(", edExpression, ",", label, ")")
        }
      }
      mainLoopHistoryEntry[["headQPR_EDExpression"]] <- edExpression

      subtree <- PCMTreeEvalNestedEDxOnTree(edExpression, tree)
      labelsSubtree <- PCMTreeGetLabels(subtree)

      if(verbose) {
        cat("Step 2.3: numTips in subtree = ", PCMTreeNumTips(subtree), "\n")
      }
      # labels of all root nodes of clades in subtree starting from its root
      cladeRootsSubtree <- partitionRootLabel

      # 2.4. listCladePartitionsSubtree: Create a list of all possible clade-partitions
      # of subtree into clades not smaller than minCladeSizes[partitionRootLevel].
      listCladePartitionsSubtree <- list()
      numPartNodes <- 1
      minCladeSizeLevel <-
        if(maxNumNodesPerCladePartition == 1) {
          # we only take one cutting node at each clade-partitioning step, so
          # no limit on the clad-size above MIN_CLADE_SIZE.
          MIN_CLADE_SIZE
        } else if(is.na(minCladeSizes[partitionRootLevel])) {
          as.integer(max(MIN_CLADE_SIZE, PCMTreeNumTips(subtree) / 8))
        } else {
          as.integer(minCladeSizes[partitionRootLevel])
        }
      while(numPartNodes <= maxNumNodesPerCladePartition) {
        listNew <- PCMTreeListCladePartitions(
          tree = subtree,
          nNodes = numPartNodes,
          minCladeSize = minCladeSizeLevel,
          tableAncestors = tableAncestors[labelsSubtree, labelsSubtree])

        if(numPartNodes == 1) {
          cladeRootsSubtree <- c(cladeRootsSubtree, labelsSubtree[unlist(listNew)])
        }

        if(length(listNew) == 0) {
          break
        } else {
          listCladePartitionsSubtree <- c(listCladePartitionsSubtree, listNew)
          numPartNodes <- numPartNodes + 1
        }
      }

      # 2.5 listCladePartitions: Unions of the nodes in bestPartition with the
      # nodes in each clade-partition of subtree; convert integer node-indices
      # in subtree to their corresponding character node-labels
      listCladePartitions <- lapply(
        listCladePartitionsSubtree,
        function(partitionSubtree) {

          partNodes <- union(bestPartition, labelsSubtree[partitionSubtree])
          PCMTreeSetRegimes(tree, as.integer(partNodes))
          partNodes2 <- PCMTreeGetStartingNodesRegimes(tree, preorder = preorderTree)

          dtTipsPerRegime <- data.table(
            regime = tree$edge.regime,
            edge2 = tree$edge[,2])[edge2 <= PCMTreeNumTips(tree),
                                   list(N=.N), keyby=regime]
          dtTipsPerRegime[, node:=nodeLabelsTree[partNodes2[regime]]]

          if(length(partNodes) > length(partNodes2) ||
             length(partNodes2) > nrow(dtTipsPerRegime) ||
             minCladeSizeLevel > dtTipsPerRegime[node%in%c(partitionRootLabel, labelsSubtree[partitionSubtree]), min(N)]) {
            NULL
          } else {
            nodeLabelsTree[partNodes2]
          }
        })


      listCladePartitions <- listCladePartitions[!sapply(listCladePartitions, is.null)]
      mainLoopHistoryEntry[["listCladePartitions"]] <- listCladePartitions

      # 2.6. Prepare "seeds" for the modelMappingIterators for each clade-partition
      if(verbose) {
        cat("Step 2.6: Preparing allowed-model-indices for the modelMappingIterators for each clade-partition\n")
      }
      PCMTreeSetRegimes(tree, as.integer(bestPartition))
      listAllowedModelTypesIndices <- list()
      for(label in union(bestPartition, cladeRootsSubtree)) {

        if(label == partitionRootLabel) {
          # we test all possible model mappings to the partition-root-table
          listAllowedModelTypesIndices[[label]] <- 1:length(modelTypes)
        } else {

          # for all other nodes in the clade-partition, we cut to the best
          # model for the clade originating at the node and the model from the
          # bestMapping over the entire tree corresponding to that node.
          bestCladeRootMapping <- tableFits[
            sapply(startingNodesRegimesLabels, function(s) match(label, s, nomatch = 0L)) == 1,
            match(mapping[[which.min(aic)]][1], modelTypes)]

          iLabel <- as.integer(label)

          # we need the if(), because PCMTreeGetRegimeForNode returns an empty vector for the root-node
          iRegime <- if(iLabel == treeRootInt) {
            1
          } else {
            PCMTreeGetRegimeForNode(tree, iLabel)
          }
          bestCladePartitionMapping <- match(bestMapping[iRegime], modelTypes)

          listAllowedModelTypesIndices[[label]] <- unique(c(bestCladePartitionMapping, bestCladeRootMapping))

          # if a previously inferred tableFits was supplied, it might contain clade-fits for model-types other than
          # the supplied modelTypes in this fit. These result in NA entries in allowedModelTypes-vectors.
          # We remove them here:
          listAllowedModelTypesIndices[[label]] <-
            listAllowedModelTypesIndices[[label]][!is.na(listAllowedModelTypesIndices[[label]])]
        }
      }

      mainLoopHistoryEntry[["listAllowedModelTypesIndices"]] <- listAllowedModelTypesIndices

      if(verbose) {
        cat("listAllowedModelTypesIndices: ",
            do.call(paste0,
                    c(as.list(capture.output(print(listAllowedModelTypesIndices))),
                      list(sep="; "))))
      }



      # 2.7 Nested foreach over the clade-partitions and the allowed model-type mappings
      if(verbose) {
        cat("level:", partitionRootLevel, "; minCladeSize:", minCladeSizeLevel,
            "; Fits on ", length(listCladePartitions), " partitions of tree into clades using the model-fits on the clades as starting points...\n")
      }

      fitsToTree <-
        foreach(
          cladePartition = listCladePartitions,
          iPartition = 1:length(listCladePartitions),
          .combine = function(...) rbindlist(list(...)),
          .multicombine = TRUE,
          .inorder = FALSE,
          .packages = (.packages()) ) %:%

        foreach(
          modelMapping = PCMIteratorMapping2(
            mapping = unlist(sapply(listAllowedModelTypesIndices[cladePartition], function(.) .[1])),
            modelTypes = seq_len(length(modelTypes)),
            allowedModelTypesIndices = lapply(listAllowedModelTypesIndices[cladePartition], sort)),

          .combine=function(...) {
            CombineTaskResults(
              ...,
              filePrefix = paste0(prefixFiles, partitionRootLevel, "_cladeparts_"),
              envNCalls = envCombineTaskResults2)
          },
          .multicombine=TRUE,
          .inorder = FALSE,
          .errorhandling = "pass",
          .packages = (.packages()) ) %op% {
            try({
              # recreate tableAncestors under a different name, to avoid exporting a potentially big
              # tableAncestors
              if(!exists("tableAncestorsTree", .GlobalEnv)) {
                if(verbose) {
                  cat("Creating tableAncestorsTree in .GlobalEnv during tree-fits\n")
                }
                assign("tableAncestorsTree", PCMTreeTableAncestors(tree, preorderTree), .GlobalEnv)
              }

              if(!exists("generatedPCMModels", .GlobalEnv) && !is.null(generatePCMModelsFun)) {
                if(verbose) {
                  cat("Calling generatePCMModelsFun()...\n")
                }
                # this should generate PCMParentClasses and PCMSpecify functions
                # for all models in modelTypes
                generatePCMModelsFun()
                assign("generatedPCMModels", TRUE, .GlobalEnv)
              }

              # don't want the names in modelMapping, the positions correspond to cladePartition
              modelMapping <- unname(modelMapping)

              # set PCMBase options from parent process: necessary if this is executed by
              # a worker process from a cluster.
              options(listPCMOptions)

              PCMTreeSetRegimes(tree, nodes = as.integer(cladePartition))

              hashCodes <- HashCodes(
                tree = tree,
                modelTypes = modelTypes,
                startingNodesRegimesLabels = cladePartition,
                modelMapping = modelMapping)
              fit <- LookupFit(tableFits = tableFits, hashCodes = hashCodes)

              if(nrow(fit) == 1 && skipFitWhenFoundInTableFits) {
                dt.row <- fit
                if(verbose) {
                  cat("Found a fit in tableFits on clade-partition: (", toString(cladePartition),
                      "); mapping: (", toString(modelMapping), ")\n")
                }
                dt.row[, duplicated:=TRUE]
              } else {
                if(verbose) {
                  cat("Performing ML fit on clade-partition: (", toString(cladePartition),
                      "); mapping: (", toString(modelMapping), ")\n")
                }

                model <- ComposeMixedGaussianFromFits(
                  tree = tree,
                  startingNodesRegimes = as.integer(cladePartition),
                  modelTypes = modelTypes,
                  k = nrow(X),
                  R = length(cladePartition),
                  mapping = modelMapping,
                  argsMixedGaussian = argsMixedGaussian,
                  tableFits = tableFits,
                  modelTypesInTableFits = modelTypesInTableFits,
                  tableAncestors = tableAncestorsTree,
                  verbose = verboseComposeMixedGaussianFromFits)

                fit <- PCMFit(
                  X = X, tree = tree, model = model, SE = SE,
                  metaI = metaIFun, positiveValueGuard = positiveValueGuard,
                  lik = lik, prior = prior, input.data = input.data,
                  config = config,
                  argsPCMParamLowerLimit = argsPCMParamLowerLimit,
                  argsPCMParamUpperLimit = argsPCMParamUpperLimit,
                  argsPCMParamLoadOrStore = argsPCMParamLoadOrStore,
                  argsConfigOptimAndMCMC = AdaptArgsConfigOptimAndMCMC(
                    model,
                    argsPCMParamLowerLimit = argsPCMParamLowerLimit,
                    argsPCMParamUpperLimit = argsPCMParamUpperLimit,
                    argsPCMParamLoadOrStore = argsPCMParamLoadOrStore,
                    argsConfigOptimAndMCMC = argsConfigOptimAndMCMC2,
                    numJitterRootRegimeFit = numJitterRootRegimeFit,
                    sdJitterRootRegimeFit = sdJitterRootRegimeFit,
                    numJitterAllRegimeFits = numJitterAllRegimeFits,
                    sdJitterAllRegimeFits = sdJitterAllRegimeFits,
                    verbose = verboseAdaptArgsConfigOptimAndMCMC),
                  verbose = verbosePCMFit)

                ll <- unname(logLik(fit))
                v_aic = unname(AIC(fit))
                vec <- c(coef(fit), logLik = ll, df = attr(ll, "df"), nobs = attr(ll, "nobs"), aic = v_aic)

                dt.row <- data.table(
                  treeEDExpression = treeEDExpression,
                  hashCodeTree = hashCodes$hashCodeTree,
                  hashCodeStartingNodesRegimesLabels = hashCodes$hashCodeStartingNodesRegimesLabels,
                  hashCodeMapping = hashCodes$hashCodeMapping,
                  startingNodesRegimesLabels = list(cladePartition),
                  mapping = list(MatchModelMapping(modelMapping, modelTypes)),
                  fitVector = list(unname(vec)),
                  logLik = ll,
                  df = attr(ll, "df"),
                  nobs = attr(ll, "nobs"),
                  aic = v_aic,
                  duplicated = FALSE)
              }

              if(printFitVectorsToConsole) {
                cat(dt.row$treeEDExpression[[1]], ":", dt.row$hashCodeTree[[1]], ":", dt.row$hashCodeStartingNodesRegimes[[1]], ":", dt.row$hashCodeMapping[[1]], ":", toString(unname(dt.row$fitVector[[1]])), "\n", sep="")
              }

              dt.row
            }, silent=FALSE)

          } # end of nested foreach body

      if(nrow(tableFits) > 0) {
        tableFits <- rbindlist(list(tableFits, fitsToTree))
      } else {
        tableFits <- fitsToTree
      }
      tableFits <- tableFits[duplicated == FALSE]
      setkey(tableFits, hashCodeTree,hashCodeStartingNodesRegimesLabels,hashCodeMapping)

      # add new entry to the main loop history
      mainLoopHistory[[length(mainLoopHistory) + 1]] <- mainLoopHistoryEntry

      # Save the current results to a file "CurrentResults<filePrefix>.RData"
      SaveCurrentResults(tableFits, mainLoopHistory, queuePartitionRoots, filePrefix = prefixFiles)

      CleanTemporaryFitFiles(filePrefix = prefixFiles)

      # 2.8 Identify the best partition and mapping after the partitioning using fitsToTree

      if( length(listCladePartitions) > 0 ) {
        if(fitsToTree[, min(aic)] < bestAIC &&
           partitionRootLevel < maxCladePartitionLevel ) {
          # a partitioning with a better AIC was found

          # put new partitioning nodes and the partition root-node in the queue with an augmented level
          queuePartitionRootsNew <- fitsToTree[, {
            iMinAIC <- which.min(aic)
            list(
              level = partitionRootLevel + 1L,
              node = as.integer(c(partitionRootLabel, setdiff(startingNodesRegimesLabels[[iMinAIC]], bestPartition))),
              partitionParentNode = as.integer(partitionRootLabel),
              hashCodeBestPartitionInitial = hashCodeStartingNodesRegimesLabels[iMinAIC],
              hashCodeBestMappingInitial = hashCodeMapping[iMinAIC],
              hashCodeBestPartitionLevel = hashCodeStartingNodesRegimesLabels[iMinAIC],
              hashCodeBestMappingLevel = hashCodeMapping[iMinAIC]
            )
          }]

          # update the best partition and mapping for the current and all upcoming entries in the queue having
          # hashCodeBestPartitionLevel == hashCodeBestPartition with the current level (partitionRootLevel)
          # doing it in this way should facilitate future inclusion of sub-optimal partitions in the queue
          hashCodeBestPartitionLevelNew <- queuePartitionRootsNew[, hashCodeBestPartitionInitial[1]]
          hashCodeBestMappingLevelNew <- queuePartitionRootsNew[, hashCodeBestMappingInitial[1]]
          queuePartitionRoots[
            c(rep(FALSE, headQPR - 1), rep(TRUE, .N - headQPR + 1)) &
              hashCodeBestPartitionLevel == hashCodeBestPartition,
            (c("hashCodeBestPartitionLevel", "hashCodeBestMappingLevel")) := list(hashCodeBestPartitionLevelNew, hashCodeBestMappingLevelNew)]
          # previously the update was done only for the nodes at the current level:
          #queuePartitionRoots[list(partitionRootLevel),
          #                    hashCodeBestPartitionLevel:=hashCodeBestPartitionLevelNew]
          #queuePartitionRoots[list(partitionRootLevel),
          #                    hashCodeBestMappingLevel:=hashCodeBestMappingLevelNew]

          # append the new nodes to the queue and reset the key
          queuePartitionRoots <- rbindlist(list(queuePartitionRoots,
                                                queuePartitionRootsNew))
          #setkey(queuePartitionRoots, level, node)
        }
      }

    } # end of main loop: while(headQPR <= nrow(queuePartitionRoots))


  } else {
    stop("ERR:04131:PCMFit:PCMFitModelMappings.R:PCMFitModelMappings:: the tree is has fewer tips than than the min clade-size in minCladeSizes (", MIN_CLADE_SIZE, "). Try with smaller minCladeSizes.")
  }

  res <- list(
    options = PCMOptions(),
    treeOriginal = treeOriginal,
    tree = tree,
    hashCodeTree = hashCodeEntireTree,
    X = X,
    SE = SE,
    tableFits = tableFits,
    queuePartitionRoots = queuePartitionRoots,
    mainLoopHistory = mainLoopHistory,
    arguments = list(
      modelTypes = modelTypes,
      metaIFun = metaIFun, generatePCMModelsFun = generatePCMModelsFun,

      positiveValueGuard = positiveValueGuard,

      lik = lik, prior = prior, input.data = input.data, config = config,

      skipFitWhenFoundInTableFits = skipFitWhenFoundInTableFits, prefixFiles = prefixFiles,

      maxCladePartitionLevel = maxCladePartitionLevel, minCladeSizes = minCladeSizes,

      listPCMOptions = listPCMOptions,

      argsMixedGaussian = argsMixedGaussian,
      argsPCMParamLowerLimit = argsPCMParamLowerLimit,
      argsPCMParamUpperLimit = argsPCMParamUpperLimit,
      argsPCMParamLoadOrStore = argsPCMParamLoadOrStore,

      argsConfigOptimAndMCMC1 = argsConfigOptimAndMCMC1,
      argsConfigOptimAndMCMC2 = argsConfigOptimAndMCMC2,

      numJitterRootRegimeFit = numJitterRootRegimeFit,
      sdJitterRootRegimeFit = sdJitterRootRegimeFit,

      numJitterAllRegimeFits = numJitterAllRegimeFits,
      sdJitterAllRegimeFits = sdJitterAllRegimeFits,

      printFitVectorsToConsole = printFitVectorsToConsole,
      setAttributes = setAttributes,

      doParallel = doParallel,

      verbose = verbose,
      verbosePCMFit = verbosePCMFit,
      verboseComposeMixedGaussianFromFits = verboseComposeMixedGaussianFromFits,
      verboseAdaptArgsConfigOptimAndMCMC = verboseAdaptArgsConfigOptimAndMCMC
    )
  )
  class(res) <- "PCMFitModelMappings"
  res
}

#' Retrieve the ML fits from the fitVectors column in a table of fits.
#' @param fitMappings an object of S3 class PCMFitModelMappings.
#' @param tableFits a data.table
#' @param modelTypesNew NULL or a character vector containing all model-types in
#'  fitMappings$arguments$modelTypes and, eventually, additional model-types.
#' @param setAttributes logical indicating if an X and tree attribute should be set to
#'   each model-object. This is used for later evaluation of the log-likelihood of the
#'   AIC coefficient for the model on the given tree and data. Using a global tree for
#'   that is a bad idea, because the model may be fit for a subtree, i.e. clade.
#'   Default FALSE.
#' @return a copy of tableFits with added column "model" and, if necessary, updated
#'  integer model-type indices in the "fitVector" column.
#' @importFrom PCMBase MixedGaussian PCMTreeEvalNestedEDxOnTree
#'
#' @export
RetrieveFittedModelsFromFitVectors <- function(
  fitMappings, tableFits = fitMappings$tableFits,
  modelTypesNew = NULL,
  setAttributes = FALSE) {

  tableFits2 <- copy(tableFits)
  tableFits2[, fittedModel:=lapply(1:.N, function(i, numRows) {
    model <- do.call(
      PCMLoadMixedGaussianFromFitVector,
      c(list(fitVector = fitVector[[i]],
           modelTypes = fitMappings$arguments$modelTypes,
           k = nrow(fitMappings$X)),
        fitMappings$arguments$argsMixedGaussian)
    )
    if(!is.null(modelTypesNew)) {
      # update the modelTypes and mapping attribute of the model appropriately
      if(is.character(modelTypesNew) && all(fitMappings$arguments$modelTypes %in% modelTypesNew) ) {
        # note that the constructor MixedGaussian accepts character vector as well as integer vector for mapping.
        mappingNew <- attr(model, "modelTypes")[attr(model, "mapping")]
        modelNew <- do.call(
          MixedGaussian,
          c(list(k = nrow(fitMappings$X),
                 modelTypes = modelTypesNew,
                 mapping = mappingNew),
            fitMappings$arguments$argsMixedGaussian))
        PCMParamLoadOrStore(modelNew, PCMParamGetShortVector(model), offset = 0, load = TRUE)
        model <- modelNew
      } else {
        stop(paste0("ERR:04141:PCMFit:PCMFitModelMappings.R:RetrieveFittedModels:: if modelTypesNew is not NULL fitMappings$arguments$modelTypes (", toString(fitMappings$arguments$modelTypes), ") should be a subset of modelTypesNew (", toString(modelTypesNew), ")."))
      }
    }
    if(setAttributes) {
      tree <- PCMTreeEvalNestedEDxOnTree(treeEDExpression[[i]], fitMappings$tree)
      PCMTreeSetRegimes(tree, match(startingNodesRegimesLabels[[i]], PCMTreeGetLabels(tree)))
      X <- fitMappings$X[, tree$tip.label]
      attr(model, "tree") <- tree
      attr(model, "X") <- X
    }
    if(numRows == 1) {
      list(model)
    } else {
      model
    }
  }, numRows = .N)]

  if(!is.null(modelTypesNew)) {
    # update the fitVectors according to the new modelTypes
    if(is.character(modelTypesNew) && all(fitMappings$arguments$modelTypes %in% modelTypesNew) ) {
      # note that the constructor MixedGaussian accepts character vector as well as integer vector for mapping.
      tableFits2[, fitVector:=lapply(1:.N, function(i) {
        treei <- PCMTreeEvalNestedEDxOnTree(treeEDExpression[[i]], fitMappings$tree)
        PCMTreeSetRegimes(treei, match(startingNodesRegimesLabels[[i]], PCMTreeGetLabels(treei)))
        par <- c(PCMGetVecParamsRegimesAndModels(fittedModel[[i]], treei), numParam = PCMParamCount(fittedModel[[i]]))
        fitVec <- fitVector[[i]]
        fitVec[1:length(par)] <- par
        fitVec
      })]
    } else {
      stop(paste0("ERR:04142:PCMFit:PCMFitModelMappings.R:RetrieveFittedModels:: if modelTypesNew is not NULL fitMappings$arguments$modelTypes (", toString(fitMappings$arguments$modelTypes), ") should be a subset of modelTypesNew (", toString(modelTypesNew), ")."))
    }
  }
  tableFits2
}

#' @export
RetrieveBestFitAIC <- function(fitMappings, rank = 1) {
  tableFits <- RetrieveFittedModelsFromFitVectors(
    fitMappings,
    fitMappings$tableFits[treeEDExpression=="tree"][order(aic)][rank],
    setAttributes = TRUE)


  res <- list(
    tree = fitMappings$tree,
    X = fitMappings$X,
    modelTypes = fitMappings$arguments$modelTypes,
    inferredRegimeNodes = tableFits$startingNodesRegimesLabels[[1]],
    inferredMapping = tableFits$mapping[[1]],
    inferredMappingIdx = match(tableFits$mapping[[1]], fitMappings$arguments$modelTypes),
    inferredModel = tableFits$fittedModel[[1]]
  )

  PCMTreeSetLabels(res$tree)
  PCMTreeSetRegimes(res$tree, res$inferredRegimeNodes)

  res[["inferredMappedModels"]] <- attr(res$inferredModel, "mapping")[res$tree$edge.regime]
  res
}

#' @importFrom data.table data.table
#' @export
PlotSearchHistory <- function(
  fit,
  sizeGreyNodepoints = 2.2, sizeColorNodepoints = 2.2,
  sizeBlackAllowedModelTypes = 1.4, sizeColorAllowedModelTypes = 1.4, sizeRankInQueue = 1.4,
  vjustBlackAllowedModelTypes = -1.6, vjustColorAllowedModelTypes = -1.6,
  ...) {
  tree <- fit$tree

  treeRootInt <- PCMTreeNumTips(tree) + 1L
  PCMTreeSetLabels(tree)

  plotList <- lapply(1:length(fit$mainLoopHistory), function(i) {
    #cat("Step ", i, "\n")

    historyEntry <- fit$mainLoopHistory[[i]]
    rootNodei <- fit$queuePartitionRoots[i, node]

    PCMTreeSetRegimes(tree, historyEntry$headQPR_Partition)

    if(length(historyEntry$listCladePartitions) > 0) {
      dtCladePartition <- data.table(node=unique(unlist(historyEntry$listCladePartitions)))
      #print(dtCladePartition[, node])
      setkey(dtCladePartition, node)

      dtCladePartition[node%in%historyEntry$headQPR_Partition, selected:=TRUE]
      dtCladePartition[!(node%in%historyEntry$headQPR_Partition), candidate:=TRUE]

      # trying to reconstruct the remaining queue is hard - better to save it at runtime;
      remainingQueuei <- fit$queuePartitionRoots[i:historyEntry$lengthQPR, list(node = as.character(node))]

      remainingQueuei[, rankInQueue:=(i + .I - 1)]
      remainingQueuei <- remainingQueuei[!is.na(node)]
      setkey(remainingQueuei, node)


      dtCladePartition[, allowedModelTypes:=sapply(node, function(n) {
          iLabel <- as.integer(n)
          # we need the if(), because PCMTreeGetRegimeForNode returns an empty vector for the root-node
          iRegime <- if(iLabel == treeRootInt) {
            historyEntry$headQPR_MappingIdx[1]
          } else {
            historyEntry$headQPR_MappingIdx[PCMTreeGetRegimeForNode(tree, iLabel)]
          }

          text <- do.call(paste,
                          c(as.list(LETTERS[unique(c(iRegime, historyEntry$listAllowedModelTypesIndices[[n]]))]),
                            list(sep="")))
          paste0("{",text,"}")
      })]
      dtCladePartition[node%in%historyEntry$headQPR_Partition, allowedModelTypesSelected:=allowedModelTypes]
      dtCladePartition[!(node%in%historyEntry$headQPR_Partition), allowedModelTypesCandidate:=allowedModelTypes]

      fitTablei <- RetrieveFittedModelsFromFitVectors(
        fit,
        LookupFit(tree,
                  fit$arguments$modelTypes,
                  historyEntry$headQPR_Mapping,
                  fit$tableFits), setAttributes = TRUE)

      ploti <- PCMTreePlot(tree, ...) %<+% as.data.frame(dtCladePartition) %<+% as.data.frame(remainingQueuei) +
        geom_nodepoint(aes(shape=selected), size=sizeColorNodepoints, na.rm = TRUE) +
        geom_nodepoint(aes(shape=candidate), size=sizeGreyNodepoints, color = "grey", na.rm = TRUE) +
        geom_text(aes(label=allowedModelTypesSelected), size=sizeColorAllowedModelTypes, vjust=vjustColorAllowedModelTypes) +
        geom_text(aes(label=allowedModelTypesCandidate), color = "black", size=sizeBlackAllowedModelTypes, vjust=vjustBlackAllowedModelTypes) +
        geom_text(aes(label=rankInQueue), color="black", size=sizeRankInQueue) +
        ggtitle(paste0("(",i,") AIC=", round(fitTablei$aic[[1]]), ", logLik=", round(fitTablei$logLik[[1]]), ", p=", fitTablei$df[[1]]))

      ploti
    } else {
      NULL
    }
  })
  plotList
}

#' @export
PlotTreeRegimesAndMapping <- function(tree, regimeNodes, mappingIdx = NULL) {
  PCMTreeSetLabels(tree)
  PCMTreeSetRegimes(tree, regimeNodes)
  pl <- PCMTreePlot(tree) #+ geom_nodelab(size = 2)

  if(!is.null(mappingIdx)) {
    dtModelMapping <- data.table(node = regimeNodes, mappedModel = LETTERS[mappingIdx])
    pl <- pl %<+% dtModelMapping +
      geom_text(aes(label = mappedModel), size=1.6, vjust=-1)
  }
  pl
}

