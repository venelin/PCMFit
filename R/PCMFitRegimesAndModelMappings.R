library(PCMBase)
library(digest)

#' Match a model mapping vector against a vector of model types
#' @param modelMapping a character or integer vector
#' @param modelTypes a vector
#' @return an integer vector or stops with an error
MatchModelMapping <- function(modelMapping, modelTypes) {
  if(is.character(modelMapping)) {
    m <- match(modelMapping, modelTypes)
    if( any(is.na(m)) ) {
      stop("ERR:04101:PCMFit:PCMFit.R:MatchModelMapping:: some of the models in modelMapping could not be matched against model-types in tableFits (", toString(modelMapping[which(is.na(m))]), ")")
    }
  } else if(is.integer(modelMapping)) {
    m <- match(modelMapping, 1:length(modelTypes))
    if( any(is.na(m)) ) {
      stop("ERR:04102:PCMFit:PCMFit.R:MatchModelMapping:: some of the integer models in modelMapping could not be matched against model-indices in tableFits (", toString(modelMapping[which(is.na(m))]), ")")
    }
  } else {
    stop("ERR:04103:PCMFit:PCMFit.R:MatchModelMapping:: modelMapping should be character or integer (",
         toString(modelMapping), ")")
  }
  m
}

#' @importFrom PCMBase PCMTreeGetStartingNodesRegimes PCMTreeGetLabels
HashCodes <- function(tree, modelTypes, modelMapping) {
  list(
    hashCodeTree = digest(PCMTreeToString(tree), serialize = FALSE),
    hashCodeStartingNodesRegimesLabels = digest(toString(PCMTreeGetLabels(tree)[PCMTreeGetStartingNodesRegimes(tree)]), serialize = FALSE),
    hashCodeMapping = digest(toString(MatchModelMapping(modelMapping, modelTypes)), serialize = FALSE)
  )
}

#' Lookup a fit vector for a given tree and model mapping in a data.table of
#' previously run fits.
#'
#' @description This is an internal function
#' @param tree a phylo object
#' @param modelTypes character vector
#' @param modelMapping an integer or character vector to by matched against modelTypes
#' @param tableFits a data.table having at least the following columns:
#' \itemize{
#' \item{hashCodeTree}{a key column of type character-vector}
#' \item{hashCodeStartingNodesRegimes}{a key column of type character-vector representing the hash-code of
#' \code{PCMTreeGetLabels(tree)[PCMTreeGetStartingNodesRegimes(tree)]}.}
#' \item{hashCodeMapping}{a key column of type character-vector}}
#' @return the corresponding fit-vector to the given tree and model mapping or
#' if no such entry is found, issues an error.
#' @importFrom digest digest
LookupFit <- function(tree, modelTypes, modelMapping, tableFits,
                      hashCodes = HashCodes(tree, modelTypes, modelMapping)) {
  tableFits[hashCodes, , mult="first", nomatch=0]
}

#' Compose an MRG model
ComposeMRGFromFits <- function(tree, startingNodesRegimes, modelTypes, k, R, mapping, argsMRG, tableFits, tableAncestors = NULL, verbose = FALSE) {
  # create an mrg model
  model <- do.call(MRG, c(list(k = k, modelTypes = modelTypes, mapping = mapping), argsMRG))

  # load listParInitOptim with the parameter vector from the ML fits to all clades
  # starting at startingNodesRegimes
  specParams <- attr(model, "specParams", exact = TRUE)
  subModelsLoaded <- rep(FALSE, length(mapping))

  for(r in 1:length(startingNodesRegimes)) {
    nr <- startingNodesRegimes[r]
    mr <- modelTypes[mapping[r]]
    tree_nr <- PCMTreeExtractClade(tree, nr, tableAncestors = tableAncestors)
    PCMTreeSetDefaultRegime(tree_nr, 1)
    fit <- LookupFit(tree = tree_nr, modelTypes = modelTypes, modelMapping = mr, tableFits = tableFits)

    if(nrow(fit) == 0) {
      stop("ERR:04120:PCMFit:PCMFit.R:ComposeMRGFromCladeFits:: no entry in tableFits for the given tree and modelMapping.")
    }

    fitVector <- fit$fitVector[[1]]

    if(verbose) {
      cat("Loading cladeFit for nr=", nr,"; mr=", mr, "\n")
    }

    fitModel <-
      do.call(PCMLoadMRGFromFitVector,
              c(list(fitVector = fitVector,
                     modelTypes = modelTypes,
                     k = k),
                argsMRG))

    for(name in names(fitModel)) {
      if(verbose) {
        cat("Setting model or global parameter with name=", name, "\n")
      }

      if(name %in% names(specParams) &&
         specParams[[name]]$type[1] %in% c("gscalar", "gvector", "gmatrix")) {
        if( specParams[[name]]$type[2] != "fixed" ) {
          # set all global non-fixed parameters to the mean of their best fits for the
          # clades
          model[[name]] <- model[[name]] + (fitModel[[name]]/R)
        }
      } else if(specParams[[name]]$type[1] == "model" &&
                class(fitModel[[name]])[1] == mr) {
        if(!subModelsLoaded[r]) {
          model[[as.character(r)]] <- fitModel[[name]]
          subModelsLoaded[r] <- TRUE
        } else {
          stop("ERR:04125:PCMFit:PCMFit.R:ComposeMRGFromCladeFits:: submodel for regime ", r, " was already loaded from a best clade fit.")
        }
      } else {
        stop("ERR:04126:PCMFit:PCMFit.R:ComposeMRGFromCladeFits:: Found a member (", name, ") in fitModel starting from node (", nr, ") and with class '", class(fitModel[[name]]), "' which is neither a global parameter nor a model of the needed type (", mr, ").")
      }
    }
  }
  model
}

AdaptArgsConfigOptimAndMCMC <- function(
  model,

  argsPCMLowerBound,
  argsPCMUpperBound,

  argsPCMSetOrGetVecParams,
  argsConfigOptimAndMCMC,
  numJitterRootRegimeFit, sdJitterRootRegimeFit,
  numJitterAllRegimeFits, sdJitterAllRegimeFits,
  verbose = FALSE
) {

  matParamsFromTableFits <- matrix(PCMGetVecParams(model), 1, PCMNumParams(model), byrow = TRUE)
  matParamsJitterRootCladeFit <- matParamsJitterAllCladeFits <- NULL

  # if there is more than one clade in the tree and numJitterRootRegimeFit > 0
  if( !is.null(model[["2"]]) && numJitterRootRegimeFit > 0 ) {
    vecParamIndex <- 1:ncol(matParamsFromTableFits)
    modelIndexParams <- model
    PCMSetOrGetVecParams(modelIndexParams, vecParamIndex)
    vecParamIndexRootClade <- as.integer(PCMGetVecParams(modelIndexParams[["1"]]))
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
      if(nrow(matParamsJitterRootCladeFit) > numJitterRootRegimeFit) {
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
      if(nrow(matParamsJitterAllCladeFits) > numJitterAllRegimeFits) {
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
      cat("Prepending the following parameter vectors to argsConfigOptimAndMCMC[['listParInitOptim']] : ")
      print(matParamsFromTableFits)
    }
    argsConfigOptimAndMCMC[["listParInitOptim"]] <-
      c(lapply(1:nrow(matParamsFromTableFits), function(i) matParamsFromTableFits[i, ]),
        argsConfigOptimAndMCMC[["listParInitOptim"]])
  } else {
    if(verbose) {
      cat("Setting listParInitOptim to : ")
      print(matParamsFromTableFits)
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

#' @importFrom foreach foreach when %do% %dopar% %:%
#' @importFrom data.table data.table rbindlist funion
#' @importFrom PCMBase PCMTreeSetLabels PCMTreeSetDefaultRegime PCMTreeEvalNestedED PCMTreeNumTips PCMTreeListCladePartitions PCMTreeToString MRG PCMOptions
#' @export
PCMFitRegimesAndModelMappingsRecursive <- function(
  X, tree, treeEDExpression = "tree", modelTypes = c("BM__noX0__noSigmae_x", "OU__noX0__noSigmae_x"),
  metaIFun = PCMInfo, positiveValueGuard = Inf,

  lik = NULL, prior = NULL, input.data = NULL, config = NULL,

  tableFits = NULL, skipFitWhenFoundInTableFits = TRUE, prefixFiles = "fits_",

  maxRecDepth = 1, minCladeSize = 25, currentRecDepth = 1,

  listPCMOptions = PCMOptions(),

  argsMRG = NULL,
  argsPCMLowerBound = NULL,
  argsPCMUpperBound = NULL,
  argsPCMSetOrGetVecParams = NULL,

  argsConfigOptimAndMCMC1 = NULL,
  argsConfigOptimAndMCMC2 = NULL,

  numJitterRootRegimeFit = 100, sdJitterRootRegimeFit = 0.5,
  numJitterAllRegimeFits = 100, sdJitterAllRegimeFits = 0.5,

  printFitVectorsToConsole = FALSE,
  setAttributes = TRUE,

  doParallel = FALSE,

  verbose = TRUE,
  verbosePCMFit = FALSE,
  verboseComposeMRGFromFits = FALSE,
  verboseAdaptArgsConfigOptimAndMCMC = FALSE
) {

  if(currentRecDepth == 1) {
    PCMTreeSetLabels(tree)
    PCMTreeSetDefaultRegime(tree, 1)
    colnames(X) <- as.character(1:PCMTreeNumTips(tree))
  } else {
    tree <- PCMTreeEvalNestedED(treeEDExpression)
    X <- X[, tree$tip.label]
  }

  if(!is.data.table(tableFits)) {
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
                           AIC = double(),
                           duplicated = logical())
  }
  setkey(tableFits, hashCodeTree,hashCodeStartingNodesRegimesLabels,hashCodeMapping )

  if(PCMTreeNumTips(tree) > minCladeSize[currentRecDepth]) {
    # start from a list containing the trivial partition into one clade equal to the whole tree
    listCladePartitions <- list(integer(0))
    numPartNodes <- 1
    while(TRUE) {
      listNew <- PCMTreeListCladePartitions(tree, numPartNodes, minCladeSize[currentRecDepth])
      if(length(listNew) == 0) {
        break
      } else {
        listCladePartitions <- c(listCladePartitions, listNew)
        numPartNodes <- numPartNodes + 1
      }
    }

    tableAncestors <- PCMTreeTableAncestors(tree)
    nodeLabels <- PCMTreeGetLabels(tree)

    `%op%` <- if(doParallel) `%dopar%` else `%do%`
    cladeRoots <- c(PCMTreeNumTips(tree) + 1,
                    unlist(listCladePartitions[sapply(listCladePartitions, length) == 1]))
    if(verbose) {
      cat("Rec-depth:", currentRecDepth, "; minCladeSize:", minCladeSize[currentRecDepth],
          "; Step 1: perform a fit on", length(cladeRoots), " clades x", length(modelTypes), " model-types...\n")
    }

    envCombineAndSave1 <- new.env()
    envCombineAndSave1$ncalls <- 0
    CombineAndSave <- function(..., filePrefix, envNCalls) {
      envNCalls$ncalls <- envNCalls$ncalls + 1
      data <- rbindlist(list(...))
      save(data, file = paste0(filePrefix, envNCalls$ncalls, ".RData"))
      data
    }

    fitsToClades <- foreach(clRoot = cladeRoots,
                            clEDExpression = paste0("E(", treeEDExpression, ",", nodeLabels[cladeRoots], ")"),
                            .combine=function(...) rbindlist(list(...)),
                            .multicombine=TRUE,
                            .inorder = FALSE,
                            .packages = (.packages()) ) %:%

      foreach(modelMapping = 1:length(modelTypes),
              .combine=function(...) {
                CombineAndSave(
                  ...,
                  filePrefix = paste0(prefixFiles, currentRecDepth, "_clades_"),
                  envNCalls = envCombineAndSave1)
              },
              .multicombine = TRUE,
              .inorder = FALSE,
              .packages = (.packages()) ) %op% {

                options(listPCMOptions)
                treeSplit <- PCMTreeSplitAtNode(tree, clRoot, tableAncestors, X)
                clade <- treeSplit$clade
                Xclade <- treeSplit$Xclade

                PCMTreeSetDefaultRegime(clade, 1)

                hashCodes <- HashCodes(tree = clade, modelTypes = modelTypes, modelMapping = modelMapping)
                fit <- LookupFit(tableFits = tableFits, hashCodes = hashCodes)

                if(nrow(fit) == 1 && skipFitWhenFoundInTableFits) {
                  dt.row <- fit
                  if(verbose) {
                    cat("Found a fit in tableFits for ", modelTypes[modelMapping], "on clade starting at node ", clRoot, " ...\n")
                  }
                  if(printFitVectorsToConsole) {
                    cat(dt.row$treeEDExpression[[1]], ":", dt.row$hashCodeTree[[1]], ":", dt.row$hashCodeStartingNodesRegimes[[1]], ":", dt.row$hashCodeMapping[[1]], ":", toString(unname(dt.row$fitVector[[1]])), "\n", sep="")
                  }
                  dt.row[, duplicated:=TRUE]

                } else {
                  if(verbose) {
                    cat("Fitting", modelTypes[modelMapping], "to clade starting at node ", clRoot, " ...\n")
                  }

                  # create an mrg model
                  model <- do.call(MRG, c(list(k = nrow(Xclade),
                                               modelTypes = modelTypes,
                                               mapping = modelMapping),
                                          argsMRG))

                  fit <- PCMFit(X = Xclade, tree = clade, model = model,
                                metaI = metaIFun, positiveValueGuard = positiveValueGuard,

                                lik = lik, prior = prior, input.data = input.data,
                                config = config,
                                argsPCMLowerBound = argsPCMLowerBound,
                                argsPCMUpperBound = argsPCMUpperBound,
                                argsPCMSetOrGetVecParams = argsPCMSetOrGetVecParams,
                                argsConfigOptimAndMCMC = argsConfigOptimAndMCMC1,
                                verbose = verbosePCMFit)

                  ll <- logLik(fit)
                  aic = AIC(fit)
                  vec <- c(coef(fit), logLik = ll, df = attr(ll, "df"), nobs = attr(ll, "nobs"), AIC = AIC(fit))
                  dt.row <- data.table(treeEDExpression = clEDExpression,
                                       hashCodeTree = hashCodes$hashCodeTree,
                                       hashCodeStartingNodesRegimesLabels = hashCodes$hashCodeStartingNodesRegimesLabels,
                                       hashCodeMapping = hashCodes$hashCodeMapping,
                                       # character names in global tree
                                       startingNodesRegimesLabels = list(PCMTreeGetLabels(clade)[PCMTreeNumTips(clade) + 1]),
                                       mapping = list(modelMapping),
                                       fitVector = list(unname(vec)),
                                       logLik = ll,
                                       df = attr(ll, "df"),
                                       nobs = attr(ll, "nobs"),
                                       AIC = aic,
                                       duplicated = FALSE)

                }
                if(printFitVectorsToConsole) {
                  cat(dt.row$treeEDExpression[[1]], ":", dt.row$hashCodeTree[[1]], ":", dt.row$hashCodeStartingNodesRegimes[[1]], ":", dt.row$hashCodeMapping[[1]], ":", toString(unname(dt.row$fitVector[[1]])), "\n", sep="")
                }
                dt.row
              }


    if(nrow(tableFits)>0) {

      tableFits <- rbindlist(list(tableFits, fitsToClades[duplicated == FALSE]))
    } else {
      tableFits <- fitsToClades
    }

    setkey(tableFits, hashCodeTree,hashCodeStartingNodesRegimesLabels,hashCodeMapping)

    SkipMapping <- function(tree, cladePartition, modelMapping, tableFits) {
      # TODO: Implement logic for skipping all but the mappings that include the
      # best fit models to each clude in cladePartition.
      FALSE
    }

    if(verbose) {
      cat("Rec-depth:", currentRecDepth, "; minCladeSize:", minCladeSize[currentRecDepth],
          "; Step 2: Fits on ", length(listCladePartitions), " partitions of tree into clades using the model-fits on the clades as starting points...\n")
    }

    envCombineAndSave2 <- new.env()
    envCombineAndSave2$ncalls <- 0
    fitsToTree <- foreach(cladePartition = listCladePartitions, iPartition = 1:length(listCladePartitions),
                             .combine=function(...) rbindlist(list(...)),
                             .multicombine = TRUE,
                             .inorder = FALSE,
                             .packages = (.packages()) ) %:%

      foreach(modelMapping = PCMIteratorMapping(rep(1L, length(cladePartition) + 1), 1:length(modelTypes)),
              .combine=function(...) {
                CombineAndSave(
                  ...,
                  filePrefix = paste0(prefixFiles, currentRecDepth, "_cladeparts_"),
                  envNCalls = envCombineAndSave2)
              },
              .multicombine=TRUE,
              .inorder = FALSE,
              .packages = (.packages()) ) %:% when( !SkipMapping(tree, cladePartition, modelMapping, tableFits) ) %op% {

                # set PCMBase options from parent process: necessary if this is executed by
                # a worker process from a cluster.
                options(listPCMOptions)

                PCMTreeSetRegimes(tree, nodes = cladePartition)
                startingNodesRegimes <- PCMTreeGetStartingNodesRegimes(tree)

                hashCodes <- HashCodes(tree = tree, modelTypes = modelTypes, modelMapping = modelMapping)
                fit <- LookupFit(tableFits = tableFits, hashCodes = hashCodes)

                if(nrow(fit) == 1 && skipFitWhenFoundInTableFits) {
                  dt.row <- fit
                  if(verbose) {
                    cat("Found a fit in tableFits for mapping: (", toString(modelMapping), ") of model-types (", modelTypes, ")\n")
                  }
                  dt.row[, duplicated:=TRUE]
                } else {
                  if(verbose) {
                    cat("Performing ML fit on mapping: (", toString(modelMapping), ") of model-types (", modelTypes, ")\n")
                    # cat("Initial model:\n", do.call(paste, c(as.list(capture.output(print(model))), sep = "\n")))
                    # cat("names(model): (", toString(names(model)), ").\n")
                    # cat("attr(model, 'regimes'): (", toString(attr(model, "regimes")), ").\n")
                    # cat("PCMRegimes(model, tree): (", toString(PCMRegimes(model, tree)), ").\n")
                  }

                  model <- ComposeMRGFromFits(tree = tree,
                                              startingNodesRegimes = startingNodesRegimes,
                                              modelTypes = modelTypes,
                                              k = nrow(X),
                                              R = length(startingNodesRegimes),
                                              mapping = modelMapping,
                                              argsMRG = argsMRG,
                                              tableFits = tableFits,
                                              tableAncestors = tableAncestors,
                                              verbose = verboseComposeMRGFromFits)

                  fit <- PCMFit(X = X, tree = tree, model = model, metaI = metaIFun, positiveValueGuard = positiveValueGuard,
                                lik = lik, prior = prior, input.data = input.data,
                                config = config,
                                argsPCMLowerBound = argsPCMLowerBound,
                                argsPCMUpperBound = argsPCMUpperBound,
                                argsPCMSetOrGetVecParams = argsPCMSetOrGetVecParams,
                                argsConfigOptimAndMCMC = AdaptArgsConfigOptimAndMCMC(
                                  model,
                                  argsPCMLowerBound = argsPCMLowerBound,
                                  argsPCMUpperBound = argsPCMUpperBound,
                                  argsPCMSetOrGetVecParams = argsPCMSetOrGetVecParams,
                                  argsConfigOptimAndMCMC = argsConfigOptimAndMCMC2,
                                  numJitterRootRegimeFit = numJitterRootRegimeFit,
                                  sdJitterRootRegimeFit = sdJitterRootRegimeFit,
                                  numJitterAllRegimeFits = numJitterAllRegimeFits,
                                  sdJitterAllRegimeFits = sdJitterAllRegimeFits,
                                  verbose = verboseAdaptArgsConfigOptimAndMCMC),
                                verbose = verbosePCMFit)

                  ll <- logLik(fit)
                  aic = AIC(fit)
                  vec <- c(coef(fit), logLik = ll, df = attr(ll, "df"), nobs = attr(ll, "nobs"), AIC = AIC(fit))

                  dt.row <- data.table(treeEDExpression = treeEDExpression,
                                       hashCodeTree = hashCodes$hashCodeTree,
                                       hashCodeStartingNodesRegimesLabels = hashCodes$hashCodeStartingNodesRegimesLabels,
                                       hashCodeMapping = hashCodes$hashCodeMapping,
                                       startingNodesRegimesLabels = list(PCMTreeGetLabels(tree)[startingNodesRegimes]),
                                       mapping = list(modelMapping),
                                       fitVector = list(unname(vec)),
                                       logLik = ll,
                                       df = attr(ll, "df"),
                                       nobs = attr(ll, "nobs"),
                                       AIC = aic,
                                       duplicated = FALSE)
                }

                if(printFitVectorsToConsole) {
                  cat(dt.row$treeEDExpression[[1]], ":", dt.row$hashCodeTree[[1]], ":", dt.row$hashCodeStartingNodesRegimes[[1]], ":", dt.row$hashCodeMapping[[1]], ":", toString(unname(dt.row$fitVector[[1]])), "\n", sep="")
                }

                dt.row
              }

    if(nrow(tableFits) > 0) {
      tableFits <- rbindlist(list(tableFits, fitsToTree[duplicated == FALSE]))

    } else {
      tableFits <- fitsToTree
    }

    setkey(tableFits, hashCodeTree,hashCodeStartingNodesRegimesLabels,hashCodeMapping)

    # save(tableFits, file=paste0(prefixFiles, currentRecDepth, "_tree_", nodeLabels[PCMTreeNumTips(tree) + 1], ".RData"))
    tableFits
  }
}

