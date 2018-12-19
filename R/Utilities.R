#' Load a mixed-regime Gaussian model from a fit vector
#' @export
PCMLoadMixedGaussianFromFitVector <- function(
  fitVector, modelTypes, k,
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

        # set all global non-fixed parameters to the mean of their best fits
        # for the clades
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

SaveCurrentResults <- function(listResults, filePrefix) {
  status <- try(
    save(listResults, file = paste0("CurrentResults_", filePrefix, ".RData")),
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
  data
}

CleanTemporaryFitFiles <- function(filePrefix) {
  status <- try(
    {
      files <- list.files(".", pattern=paste0("^", filePrefix, ".*.RData"))
      if(length(files) > 0) {
        do.call(file.remove, as.list(files))
      }

    }, silent = TRUE)
  if(class(status) == 'try-error') {
    warning(paste0("An error occurred while cleaning up temporary data-files: ", status))
  }
}

InitTableFits <- function(
  modelTypes,
  fitMappingsPrev = NULL,
  tableFitsPrev = fitMappingsPrev$tableFits,
  modelTypesInTableFitsPrev = NULL,
  verbose = FALSE) {

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
  attr(tableFits, "modelTypes") <- modelTypesInTableFits
  tableFits
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


