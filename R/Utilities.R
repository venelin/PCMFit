#' Load a mixed-regime Gaussian model from a fit vector
#' @export
PCMLoadMixedGaussianFromFitVector <- function(
  fitVector, modelTypes, k,
  remapModelTypeIndicesInFitVector = seq_along(modelTypes), ...) {
  # the last entries in fitVector are in the following order from left to right:
  # numNumericParams, logLik, df, nobs, score;
  # the first elements from 1 to numNumericParams are the actual numeric
  # parameters; the entries that follow are between the numeric parameters and
  # numNumericParams. These must be a pair number 2R, where R is the number of
  # regimes in the tree. The first R of these are the starting nodes for the R
  # regimes (starting from the root-node, which is always the starting node for
  # the first regime); The second R are the indices that map modelTypes to the
  # regimes on the tree.
  last <- length(fitVector)
  v_score <- fitVector[last]
  nobs <- fitVector[last - 1]
  df <- fitVector[last - 2]
  ll <- fitVector[last - 3]
  numNumericParams <- fitVector[last - 4]

  mappingModelsToRegimes <- matrix(
    as.integer(fitVector[(numNumericParams+1):(last-5)]), nrow = 2, byrow = TRUE)

  model <- MixedGaussian(
    k = k, modelTypes = modelTypes,
    mapping = remapModelTypeIndicesInFitVector[mappingModelsToRegimes[2, ]], ...)

  PCMParamLoadOrStore(model, fitVector, offset = 0, load = TRUE)

  attr(model, "startingNodesRegimes") <- mappingModelsToRegimes[1,]
  attr(model, "ll") <- ll
  attr(model, "df") <- df
  attr(model, "nobs") <- nobs
  attr(model, "score") <- v_score
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
      stop("MatchModelMapping:: some of the models in modelMapping could not be matched against model-types in tableFits (", toString(modelMapping[which(is.na(m))]), ")")
    }
  } else if( is.integer(modelMapping) ) {
    m <- match(modelMapping, seq_along(modelTypes))
    if( any(is.na(m)) ) {
      stop("MatchModelMapping:: some of the integer models in modelMapping could not be matched against model-indices in tableFits (", toString(modelMapping[which(is.na(m))]), ")")
    }
  } else {
    stop(
      "MatchModelMapping:: modelMapping should be character or integer (",
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
  model <- do.call(
    MixedGaussian,
    c(list(k = k, modelTypes = modelTypes, mapping = mapping), argsMixedGaussian))

  # load listParInitOptim with the parameter vector from the ML fits to all clades
  # starting at startingNodesRegimes
  spec <- attr(model, "spec", exact = TRUE)
  subModelsLoaded <- rep(FALSE, length(mapping))

  for(r in seq_len(length(startingNodesRegimes))) {
    nr <- startingNodesRegimes[r]
    mr <- modelTypes[mapping[r]]
    tree_nr <- PCMTreeExtractClade(tree, nr, tableAncestors = tableAncestors)
    PCMTreeSetPartition(tree_nr)
    # model_mr <- do.call(
    #   MixedGaussian,
    #   c(list(k = k, modelTypes = modelTypes, mapping = mapping[r]), argsMixedGaussian))

    fit <- LookupFit(
      tree = tree_nr, modelTypes = modelTypes, modelMapping = mr,
      tableFits = tableFits)
    #fit <- LookupFit2(tree = tree_nr, model = model_nr, tableFits = tableFits)

    if(nrow(fit) == 0) {
      stop("ComposeMixedGaussianFromCladeFits:: no entry in tableFits for the given tree and modelMapping.")
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
          stop(
            paste0(
              "ComposeMixedGaussianFromCladeFits:: a global parameter ",
              name,
              " in a fit-model has a different short-vector length from the to-be fit model length(vecToFit)=",
              length(vecCurrent), ", length(vecFitModel)=", length(vecFitModel),
              "; class model to fit = ", mr,
              "; class fitted model = ", toString(class(fitModel[[name]]))))
        }
        vecCurrent <- vecCurrent + vecFitModel/R
        PCMParamLoadOrStore(model[[name]], vecCurrent, 0, k, 1, load = TRUE)

      } else if(is.PCM(spec[[name]]) && class(fitModel[[name]])[1] == mr) {
        if(!subModelsLoaded[r]) {
          model[[as.character(r)]] <- fitModel[[name]]
          subModelsLoaded[r] <- TRUE
        } else {
          stop(
            "ComposeMixedGaussianFromCladeFits:: submodel for regime ",
            r, " was already loaded from a best clade fit.")
        }
      } else {
        stop(
          "ComposeMixedGaussianFromCladeFits:: Found a member (",
          name, ") in fitModel starting from node (", nr, ") and with class '",
          class(fitModel[[name]]),
          "' which is neither a global parameter nor a model of the needed type (",
          mr, ").")
      }
    }
  }
  model
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

#' @export
SaveTempWorkerResults <- function(fitsNew, filePrefix) {
  workerPid <- Sys.getpid()
  fileName <- paste0(filePrefix, "_worker_", workerPid, ".RData")

  # load file with fits table for this worker
  status <- try({
    # previous fits stored in file
    fits <- NULL
    # loads a variable fits
    if(file.exists(fileName)) {
      load(fileName)
    }
    fits <- rbindlist(list(fits, fitsNew))
    save(fits, file = fileName)
  }, silent = TRUE)
}

# combine fits from parallel tasks into a data.table and saves this data.table to
# an .RData file
CombineTaskResults <- function(..., envNCalls) {
  envNCalls$ncalls <- envNCalls$ncalls + 1
  data <- rbindlist(
    lapply(
      list(...),
      function(dt) {
        if(is.data.table(dt)) {
          dt
        } else {
          cat("ERROR in foreach task: no data.table returned: \n", toString(dt),
              "\n")
          NULL
        }
      }),
    use.names = TRUE)
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

  if(!is.null(fitMappingsPrev)) {
    tableFits <- tableFitsPrev <- fitMappingsPrev$tableFits
    if(!identical(modelTypes, fitMappingsPrev$arguments$modelTypes)) {

      # this should remap the model-type indices in the fit vectors, show
      # table fits is correctly converted.
      tableFits <- RetrieveFittedModelsFromFitVectors(
        fitMappings = fitMappingsPrev, tableFits = tableFitsPrev,
        modelTypesNew = modelTypes)
    }
  }

  if(!is.data.table(tableFits)) {
    if(verbose) {
      cat("Initiating tableFits...\n")
    }
    tableFits = data.table(hashCodeTree = character(),
                           hashCodeStartingNodesRegimesLabels = character(),
                           hashCodeMapping = character(),
                           treeEDExpression = character(),
                           startingNodesRegimesLabels = list(),
                           mapping = list(),
                           fitVector = list(),
                           logLik = double(),
                           df = integer(),
                           nobs = integer(),
                           score = double(),
                           duplicated = logical())
  }

  setkey(tableFits, hashCodeTree,hashCodeStartingNodesRegimesLabels,hashCodeMapping )
  attr(tableFits, "modelTypes") <- modelTypes
  tableFits
}

#' @importFrom data.table is.data.table
#' @export
UpdateTableFits <- function(tableFits, newFits) {
  if(!is.data.table(newFits) && !is.data.table(tableFits)) {
    stop("Both newFits and tableFits are not data.table objects!")
  } else if(!is.data.table(newFits)) {
    # swap the two arguments
    newFits2 <- tableFits
    tableFits <- newFits
    newFits <- newFits2
  }
  if(is.null(tableFits) || !is.data.table(tableFits) || nrow(tableFits) == 0) {
    tableFits <- newFits
  } else {
    #nrow(tableFits) > 0
    tableFits <- rbindlist(list(newFits, tableFits), use.names = TRUE)
  }

  tableFits[, .SD[which.min(score)], keyby = list(
      hashCodeTree, hashCodeStartingNodesRegimesLabels, hashCodeMapping)]
}

#' Retrieve the ML fits from the fitVectors column in a table of fits.
#' @param fitMappings an object of S3 class PCMFitModelMappings.
#' @param tableFits a data.table
#' @param modelTypesNew NULL or a character vector containing all model-types in
#'  fitMappings$arguments$modelTypes and, eventually, additional model-types.
#' @param setAttributes logical indicating if an X and tree attribute should be
#' set to each model-object. This is used for later evaluation of the
#' log-likelihood of the score for the model on the given tree and
#' data. Using a global tree for that is a bad idea, because the model may be
#' fit for a subtree, i.e. clade. Default FALSE.
#' @return a copy of tableFits with added column "model" and, if necessary,
#' updated integer model-type indices in the "fitVector" column.
#' @importFrom PCMBase MixedGaussian PCMTreeEvalNestedEDxOnTree
#' @importFrom data.table setnames
#' @export
RetrieveFittedModelsFromFitVectors <- function(
  fitMappings,

  tableFits = fitMappings$tableFits,
  modelTypes = fitMappings$arguments$modelTypes,
  modelTypesNew = NULL,
  argsMixedGaussian = fitMappings$arguments$argsMixedGaussian,

  X = fitMappings$X,
  tree = fitMappings$tree,
  SE = fitMappings$SE,

  setAttributes = FALSE) {

  if(is.null(tableFits$score)) {
    # in previous versions the column score was named aic
    setnames(tableFits, "aic", "score")
  }

  # Copy all arguments into a list
  # We establish arguments$<argument-name> as a convention for accessing the
  # original argument value.
  arguments <- as.list(environment())

  tableFits2 <- copy(tableFits)
  tableFits2[, fittedModel:=lapply(seq_len(.N), function(i, numRows) {
    model <- do.call(
      PCMLoadMixedGaussianFromFitVector,
      c(list(fitVector = fitVector[[i]],
             modelTypes = modelTypes,
             k = nrow(X)),
        argsMixedGaussian)
    )
    if(!is.null(modelTypesNew)) {
      # update the modelTypes and mapping attribute of the model appropriately
      if(is.character(modelTypesNew) && all(modelTypes %in% modelTypesNew) ) {
        # note that the constructor MixedGaussian accepts character vector as well as integer vector for mapping.
        mappingNew <- attr(model, "modelTypes")[attr(model, "mapping")]
        modelNew <- do.call(
          MixedGaussian,
          c(list(k = nrow(X),
                 modelTypes = modelTypesNew,
                 mapping = mappingNew),
            argsMixedGaussian))
        PCMParamLoadOrStore(modelNew, PCMParamGetShortVector(model), offset = 0, load = TRUE)
        model <- modelNew
      } else {
        stop(
          paste0(
            "RetrieveFittedModels:: if modelTypesNew is not NULL fitMappings$arguments$modelTypes (",
            toString(modelTypes),
            ") should be a subset of modelTypesNew (",
            toString(modelTypesNew), ")."))
      }
    }
    if(setAttributes) {
      tree <- PCMTreeEvalNestedEDxOnTree(
        treeEDExpression[[i]], PCMTree(arguments$tree))
      PCMTreeSetPartition(
        tree, match(startingNodesRegimesLabels[[i]], PCMTreeGetLabels(tree)))
      X <- arguments$X[, tree$tip.label]
      SE <- arguments$SE[, tree$tip.label]
      if(is.null(SE)) {
        SE <- X
        SE[] <- 0.0
      }
      attr(model, "tree") <- tree
      attr(model, "X") <- X
      attr(model, "SE") <- SE
    }
    if(numRows == 1) {
      list(model)
    } else {
      model
    }
  }, numRows = .N)]

  if(!is.null(modelTypesNew)) {
    # update the fitVectors according to the new modelTypes
    if(is.character(modelTypesNew) && all(modelTypes %in% modelTypesNew) ) {
      # note that the constructor MixedGaussian accepts character vector as well as integer vector for mapping.
      tableFits2[, fitVector:=lapply(seq_len(.N), function(i) {
        treei <- PCMTreeEvalNestedEDxOnTree(
          treeEDExpression[[i]], PCMTree(arguments$tree))
        PCMTreeSetPartition(
          treei, match(startingNodesRegimesLabels[[i]], PCMTreeGetLabels(treei)))
        par <- c(PCMGetVecParamsRegimesAndModels(fittedModel[[i]], treei), numParam = PCMParamCount(fittedModel[[i]]))
        fitVec <- fitVector[[i]]
        fitVec[seq_len(length(par))] <- par
        fitVec
      })]
    } else {
      stop(paste0("RetrieveFittedModels:: if modelTypesNew is not NULL fitMappings$arguments$modelTypes (", toString(modelTypes), ") should be a subset of modelTypesNew (", toString(modelTypesNew), ")."))
    }
  }
  tableFits2
}


#' @importFrom data.table setnames
#' @export
RetrieveBestFitScore <- function(fitMappings, rank = 1) {
  if(is.null(fitMappings$tableFits$score)) {
    # the fit was produced with a previous version where the score column
    # was named aic.
    setnames(fitMappings$tableFits, old = "aic", new = "score")
  }
  tableFits <- RetrieveFittedModelsFromFitVectors(
    fitMappings = fitMappings,
    tableFits = fitMappings$tableFits[treeEDExpression=="tree"][order(score)][rank],
    setAttributes = TRUE)


  res <- list(
    tree = PCMTree(fitMappings$tree),
    X = fitMappings$X,
    modelTypes = fitMappings$arguments$modelTypes,
    inferredRegimeNodes = tableFits$startingNodesRegimesLabels[[1]],
    inferredMapping = tableFits$mapping[[1]],
    inferredMappingIdx = match(tableFits$mapping[[1]], fitMappings$arguments$modelTypes),
    inferredModel = tableFits$fittedModel[[1]]
  )

  PCMTreeSetLabels(res$tree)
  PCMTreeSetPartition(res$tree, res$inferredRegimeNodes)

  res[["inferredMappedModels"]] <- attr(res$inferredModel, "mapping")[res$tree$edge.regime]
  res
}

#' @importFrom PCMBase PCMDefaultObject PCMParamSetByName
#' @export
LearnCladeFitsFromSubmodels <- function(
  cladeFits,
  modelTypes,
  subModels,
  argsMixedGaussian,
  metaIFun = PCMInfo,
  scoreFun,
  X, tree, SE,
  verbose = FALSE) {

  cladeFitsNew <- cladeFits[integer(0L)]
  count <- 0L
  cladeRoots <- c()
  listAllowedModelTypesIndices <- list()

  for(modelType in names(subModels)) {
    for(edExpr in unique(cladeFits[, treeEDExpression])) {
      subModelType <- subModels[modelType]

      cladeFits2 <- cladeFits[treeEDExpression == edExpr][
        unlist(mapping) %in% modelTypes[c(modelType, subModelType)]]

      cladeFits2[, modelTypeName:=names(modelTypes)[match(unlist(mapping), modelTypes)]]
      setkey(cladeFits2, modelTypeName)

      if(nrow(cladeFits2) == 2L &&
         nrow(cladeFits2[list(modelType)]) == 1L &&
         nrow(cladeFits2[list(subModelType)]) == 1L &&
         cladeFits2[list(modelType), logLik] < cladeFits2[list(subModelType), logLik]) {
        count <- count+1
        if(verbose) {
          cat(
            count, ". ",
            'treeEDExpr=', edExpr,
            ': substituting parameters for modelType=', modelType,
            '(ll=', round(cladeFits2[list(modelType), logLik], 2), ')',
            ' with parameters from subModelType=', subModelType,
            '(ll=', round(cladeFits2[list(subModelType), logLik], 2), ')', '\n')
        }

        cladeFits2Models <- RetrieveFittedModelsFromFitVectors(
          fitMappings = NULL,
          tableFits = cladeFits2,
          modelTypes = modelTypes,
          modelTypesNew = NULL,
          argsMixedGaussian = argsMixedGaussian,

          X = X,
          tree = tree,
          SE = SE,

          setAttributes = TRUE)

        model <- cladeFits2Models[list(modelType), fittedModel[[1]]]
        subModel <- cladeFits2Models[list(subModelType), fittedModel[[1]]]
        model2 <- PCMDefaultObject(spec = attr(model, 'spec'), model = model)
        attributes(model2) <- attributes(model)
        attr(model2, "PCMInfoFun") <- metaIFun(
          X = attr(model2, "X", exact = TRUE),
          tree = attr(model2, "tree", exact = TRUE),
          model = model2,
          SE = attr(model2, "SE", exact = TRUE))
        PCMParamSetByName(model2, subModel, inplace = TRUE, deepCopySubPCMs = TRUE)

        vecModel <- cladeFits2Models[list(modelType)]$fitVector[[1]]
        idxParams <- seq_len(PCMParamCount(model))
        idxLogLik <- length(vecModel) - 3
        idxScore <- length(vecModel)
        vecModel[idxParams] <- unname(PCMParamGetShortVector(model2))
        vecModel[idxLogLik] <- unname(logLik(model2))
        vecModel[idxScore] <- unname(scoreFun(model2))

        cladeFitsNewEntry <- cladeFits2[list(modelType)]
        cladeFitsNewEntry[, modelTypeName:=NULL]

        cladeFitsNewEntry[,fitVector:=list(list(vecModel))]
        cladeFitsNewEntry[,logLik:=vecModel[[idxLogLik]]]
        cladeFitsNewEntry[,score:=vecModel[[idxScore]]]

        cladeFitsNew <- rbindlist(list(cladeFitsNew, cladeFitsNewEntry))

        cladeRoot <- cladeFitsNewEntry$startingNodesRegimesLabels[[1L]]
        cladeRoots <- c(cladeRoots, cladeRoot)
        listAllowedModelTypesIndices[[as.character(cladeRoot)]] <- c(
          listAllowedModelTypesIndices[[as.character(cladeRoot)]],
          match(modelTypes[modelType], modelTypes))
      }
    }
  }
  if(!is.null(cladeFitsNew)) {
    setkey(
      cladeFitsNew,
      hashCodeTree, hashCodeStartingNodesRegimesLabels, hashCodeMapping)
  }
  list(
    cladeFitsNew = cladeFitsNew,
    listPartitions = as.list(unique(cladeRoots)),
    listAllowedModelTypesIndices = listAllowedModelTypesIndices)
}

#' @importFrom PCMBase PCMTreeGetPartition PCMTreeGetLabels
#' @importFrom digest digest
#' @export
HashCodes <- function(
  tree, modelTypes, startingNodesRegimesLabels, modelMapping) {

  orderPNLs <- order(as.integer(startingNodesRegimesLabels))
  list(
    hashCodeTree = digest(PCMTreeToString(tree), serialize = FALSE),
    hashCodeStartingNodesRegimesLabels = digest(
      toString(startingNodesRegimesLabels[orderPNLs]), serialize = FALSE),
    hashCodeMapping = digest(
      toString(MatchModelMapping(modelMapping[orderPNLs], modelTypes)), serialize = FALSE)
  )
}

#' Lookup a fit vector for a given tree and model mapping in a data.table of
#' previously run fits.
#'
#' @param tree a phylo object
#' @param modelTypes character vector
#' @param modelMapping an integer or character vector to be matched against modelTypes
#' @param tableFits a data.table having at least the following columns:
#' \itemize{
#' \item{hashCodeTree}{an MD5 key column of type character-vector}
#' \item{hashCodePartitionNodeLabels}{an MD5 key column of type character-vector
#'  representing the hash-code of
#'  \code{PCMTreeGetLabels(tree)[PCMTreeGetPartition(tree)]}.}
#' \item{hashCodeMapping}{an MD5 key column of type character-vector}}
#' @return the corresponding fit-vector to the given tree and model mapping or
#' if no such entry is found, issues an error.
#' @importFrom digest digest
#' @importFrom data.table setnames
#' @export
LookupFit <- function(
  tree, modelTypes, modelMapping, tableFits,
  hashCodes = HashCodes(tree = tree,
                        modelTypes = modelTypes,
                        startingNodesRegimesLabels =
                          PCMTreeGetLabels(tree)[PCMTreeGetPartition(tree)],
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
  tree <- PCMTree(fit$tree)

  treeRootInt <- PCMTreeNumTips(tree) + 1L
  PCMTreeSetLabels(tree)

  plotList <- lapply(1:length(fit$mainLoopHistory), function(i) {
    #cat("Step ", i, "\n")

    historyEntry <- fit$mainLoopHistory[[i]]
    rootNodei <- fit$queuePartitionRoots[i, node]

    PCMTreeSetPartition(tree, historyEntry$headQPR_Partition)

    if(length(historyEntry$listPartitions) > 0) {
      dtCladePartition <- data.table(node=unique(unlist(historyEntry$listPartitions)))
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
        # we need the if(), because PCMTreeGetPartsForNodes returns an empty
        # vector for the root-node
        iRegime <- if(iLabel == treeRootInt) {
          historyEntry$headQPR_MappingIdx[1]
        } else {
          historyEntry$headQPR_MappingIdx[
            PCMTreeGetPartsForNodes(tree, iLabel)]
        }

        text <- do.call(
          paste,
          c(as.list(LETTERS[
            unique(c(iRegime, historyEntry$listAllowedModelTypesIndices[[n]]))]),
            list(sep="")))
        paste0("{",text,"}")
      })]

      dtCladePartition[
        node%in%historyEntry$headQPR_Partition,
        allowedModelTypesSelected:=allowedModelTypes]
      dtCladePartition[
        !(node%in%historyEntry$headQPR_Partition),
        allowedModelTypesCandidate:=allowedModelTypes]

      fitTablei <- RetrieveFittedModelsFromFitVectors(
        fit,
        LookupFit(
          tree = tree,
          modelTypes = fit$arguments$modelTypes,
          modelMapping = historyEntry$headQPR_Mapping,
          tableFits = fit$tableFits), setAttributes = TRUE)

      ploti <- PCMTreePlot(tree, ...) %<+% as.data.frame(dtCladePartition) %<+% as.data.frame(remainingQueuei) +
        geom_nodepoint(aes(shape=selected), size=sizeColorNodepoints, na.rm = TRUE) +
        geom_nodepoint(aes(shape=candidate), size=sizeGreyNodepoints, color = "grey", na.rm = TRUE) +
        geom_text(aes(label=allowedModelTypesSelected), size=sizeColorAllowedModelTypes, vjust=vjustColorAllowedModelTypes) +
        geom_text(aes(label=allowedModelTypesCandidate), color = "black", size=sizeBlackAllowedModelTypes, vjust=vjustBlackAllowedModelTypes) +
        geom_text(aes(label=rankInQueue), color="black", size=sizeRankInQueue) +
        ggtitle(paste0("(",i,") score=", round(fitTablei$score[[1]]), ", logLik=", round(fitTablei$logLik[[1]]), ", p=", fitTablei$df[[1]]))

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
  PCMTreeSetPartition(tree, regimeNodes)
  pl <- PCMTreePlot(tree) #+ geom_nodelab(size = 2)

  if(!is.null(mappingIdx)) {
    dtModelMapping <- data.table(node = regimeNodes, mappedModel = LETTERS[mappingIdx])
    pl <- pl %<+% dtModelMapping +
      geom_text(aes(label = mappedModel), size=1.6, vjust=-1)
  }
  pl
}


