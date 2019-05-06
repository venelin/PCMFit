#' @importFrom utils tail
PCMFitModelMappingsToCladePartitions <- function(
  X, tree, modelTypes,
  SE = matrix(0.0, nrow(X), PCMTreeNumTips(tree)),

  scoreFun = AIC,

  listPartitions = NULL,
  listHintModels = rep(list(NULL), length(listPartitions)),
  listNamesInHintModels = rep(list(character(0)), length(listPartitions)),

  listAllowedModelTypesIndices = NULL,

  fitClades = FALSE,

  generatePCMModelsFun = NULL,
  metaIFun = PCMInfo, positiveValueGuard = Inf,

  fitMappingsPrev = NULL,
  tableFitsPrev = fitMappingsPrev$tableFits,
  modelTypesInTableFitsPrev = NULL,

  skipFitWhenFoundInTableFits = TRUE,

  argsMixedGaussian = NULL,

  argsConfigOptim = NULL,


  listPCMOptions = PCMOptions(),
  doParallel = FALSE,

  preorderTree = NULL,
  tableAncestors = NULL,

  prefixFiles = "fits_",
  saveTempWorkerResults = TRUE,
  printFitVectorsToConsole = FALSE,

  verbose = TRUE,
  debug = FALSE
) {

  tree <- PCMTree(tree)
  treeEDExpression = "tree"

  tableFits <- InitTableFits(modelTypes,
                             fitMappingsPrev,
                             tableFitsPrev,
                             modelTypesInTableFitsPrev,
                             verbose = verbose)

  modelTypesInTableFits <- attr(tableFits, "modelTypes")

  if(is.null(preorderTree)) {
    preorderTree <- PCMTreePreorder(tree)
  }

  if(is.null(tableAncestors)) {
    tableAncestors <- PCMTreeTableAncestors(tree, preorder = preorderTree)
  }

  `%op%` <- if(isTRUE(doParallel) ||
               (is.numeric(doParallel) && doParallel > 1)) `%dopar%` else `%do%`

  envCombineTaskResults <- new.env()
  envCombineTaskResults$ncalls <- 0

  EDExpressions <- if(fitClades) {
    nodeLabelsTree <- PCMTreeGetLabels(tree)
    cladeRoots <- sapply(listPartitions, function(.) .[1])
    paste0("E(", treeEDExpression, ",", nodeLabelsTree[cladeRoots], ")")
  } else {
    rep(treeEDExpression, length(listPartitions))
  }

  GetAllowedModelTypesForPartition <- function(
    listAllowedModelTypesIndices, partition, iPartition) {

    if( is.null(names(listAllowedModelTypesIndices)) ) {
      if(length(listAllowedModelTypesIndices) < iPartition) {
        stop(
          paste(
            "GetAllowedModelTypesForPartition: ",
            "length(listAllowedModelTypesIndices) = ",
            length(listAllowedModelTypesIndices), "\n",
            "iPartition = ", iPartition, "\n",
            "iPartition is bigger than length(listAllowedModelTypesIndices);",
            "listAllowedModelTypesIndices should either be a named list",
            "with names corresponding to node-labels in the tree or",
            " an unnamed list containing one such named list for each entry",
            " in listPartitions."))
      } else {
        listAllowedModelTypesIndices[[iPartition]][as.character(partition)]
      }
    } else {
      listAllowedModelTypesIndices[as.character(partition)]
    }
  }

  fits <-
    foreach(
      cladePartition = listPartitions,
      iCladePartition = seq_along(listPartitions),
      hintModel = listHintModels,
      namesInHintModel = listNamesInHintModels,
      EDExpression = EDExpressions,
      .combine = function(...) {
        listDTs <- list(...)
        for(i in seq_along(listDTs)) {
          if(!is.data.table(listDTs[[i]])) {
            cat(
              "PCMFitModelMappingsToPartitions: Error in .combine function \n",
              "for outer foreach loop: item i=", i, " from task results of ",
              "innter foreach loop was not a data.table, but was as follows:\n")
            print(listDTs[[i]])
          }
        }
        rbindlist(listDTs, use.names = TRUE)
      },
      .multicombine = TRUE,
      .inorder = FALSE,
      .packages = (.packages()) ) %:%

    foreach(
      modelMapping = PCMIteratorMapping2(
        mapping = unlist(
          sapply(
            GetAllowedModelTypesForPartition(
              listAllowedModelTypesIndices,
              cladePartition,
              iCladePartition),
            function(.) .[1])),
        modelTypes = seq_along(modelTypes),
        allowedModelTypesIndices = lapply(
          GetAllowedModelTypesForPartition(
            listAllowedModelTypesIndices,
            cladePartition,
            iCladePartition),
          sort)),

      .combine=function(...) {
        CombineTaskResults(
          ...,
          envNCalls = envCombineTaskResults)
      },
      .multicombine = TRUE,
      .inorder = FALSE,
      .errorhandling = "pass",
      .packages = (.packages()) ) %op% {
        try({
          # recreate tableAncestors under a different name, to avoid exporting a
          # potentially big tableAncestors
          if(!exists("tableAncestorsTree", .GlobalEnv)) {
            if(verbose) {
              cat("Creating tableAncestorsTree in .GlobalEnv during tree-fits\n")
            }
            assign("tableAncestorsTree",
                   PCMTreeTableAncestors(tree, preorderTree), .GlobalEnv)
          }

          if(!exists("generatedPCMModels", .GlobalEnv) &&
             !is.null(generatePCMModelsFun)) {
            if(verbose) {
              cat("Calling generatePCMModelsFun()...\n")
            }
            # this should generate PCMParentClasses and PCMSpecify functions
            # for all models in modelTypes
            generatePCMModelsFun()
            assign("generatedPCMModels", TRUE, .GlobalEnv)
          }

          # don't want the names in modelMapping, the positions correspond to
          # cladePartition
          modelMapping <- unname(modelMapping)

          # set PCMBase options from parent process: necessary if this is
          # executed by a worker process from a cluster.
          options(listPCMOptions)

          if(fitClades) {
            k <- nrow(X)
            XSE <- rbind(X, SE)
            treeSplit <- PCMTreeSplitAtNode(
              tree, as.integer(cladePartition[1]), tableAncestorsTree, XSE)
            treeForFit <- treeSplit$clade
            XForFit <- treeSplit$Xclade[1:k, , drop = FALSE]
            SEForFit <- treeSplit$Xclade[(k+1):(2*k), , drop = FALSE]
            PCMTreeSetPartition(treeForFit)

            # create an MixedGaussian model
            modelForFit <- do.call(MixedGaussian,
                                   c(list(k = nrow(XForFit),
                                          modelTypes = modelTypes,
                                          mapping = modelMapping),
                                     argsMixedGaussian))
          } else {
            treeForFit <- tree
            PCMTreeSetPartition(treeForFit, nodes = as.integer(cladePartition))
            XForFit <- X
            SEForFit <- SE

            modelForFit <- ComposeMixedGaussianFromFits(
              tree = treeForFit,
              startingNodesRegimes = as.integer(cladePartition),
              modelTypes = modelTypes,
              k = nrow(XForFit),
              R = length(cladePartition),
              mapping = modelMapping,
              argsMixedGaussian = argsMixedGaussian,
              tableFits = tableFits,
              modelTypesInTableFits = modelTypesInTableFits,
              tableAncestors = tableAncestorsTree,
              verbose = debug)

            if(!is.null(hintModel)) {
              PCMParamSetByName(modelForFit, hintModel[namesInHintModel])
            }
          }

          hashCodes <- HashCodes(
            tree = treeForFit,
            modelTypes = modelTypes,
            startingNodesRegimesLabels = cladePartition,
            modelMapping = modelMapping)
          fit <- LookupFit(tableFits = tableFits, hashCodes = hashCodes)

          if(nrow(fit) == 1 && skipFitWhenFoundInTableFits) {
            dt.row <- fit
            if(verbose) {
              cat(
                "Found a fit in tableFits on '",
                EDExpression, "', partition: (", toString(cladePartition),
                  "), mapping: (", toString(names(modelTypes)[modelMapping]),
                ")\n", sep = "")
            }
            dt.row[, duplicated:=TRUE]
          } else {
            if(verbose) {
              cat(
                "Performing ML fit on '",
                EDExpression, "', partition: (", toString(cladePartition),
                "), mapping: (", toString(names(modelTypes)[modelMapping]),
                ")\n", sep = "")
            }

            if(fitClades) {
              matParInitRunif <- PCMParamRandomVecParams(
                o = modelForFit,
                k = PCMNumTraits(modelForFit),
                R = PCMNumRegimes(modelForFit),
                n = argsConfigOptim$numRunifInitVecParams,
                argsPCMParamLowerLimit = argsConfigOptim$argsPCMParamLowerLimit,
                argsPCMParamUpperLimit = argsConfigOptim$argsPCMParamUpperLimit)
              matParInitGuess <- GuessInitVecParams(
                o = modelForFit,
                n = 1L,
                argsPCMParamLowerLimit = argsConfigOptim$argsPCMParamLowerLimit,
                argsPCMParamUpperLimit = argsConfigOptim$argsPCMParamUpperLimit,
                X = X, tree = tree, SE = SE,
                varyParams = FALSE)
              matParInitGuessVaryParams <- GuessInitVecParams(
                o = modelForFit,
                n = argsConfigOptim$numGuessInitVecParams,
                argsPCMParamLowerLimit = argsConfigOptim$argsPCMParamLowerLimit,
                argsPCMParamUpperLimit = argsConfigOptim$argsPCMParamUpperLimit,
                X = X, tree = tree, SE = SE,
                varyParams = TRUE)

              matParInit <- rbind(matParInitRunif,
                                  matParInitGuess,
                                  matParInitGuessVaryParams)
              if(nrow(fit) == 1L) {
                cat("Completing matParInit with a vector from previous fit.\n")
                matParInit <- rbind(
                  matParInit, fit$fitVector[[1]][seq_len(ncol(matParInit))])
                print(matParInit[tail(nrow(matParInit)),])
              }

              fit <- PCMFit(
                X = XForFit, tree = treeForFit, model = modelForFit,
                SE = SEForFit,
                metaI = metaIFun, positiveValueGuard = positiveValueGuard,
                matParInit = matParInit,
                numCallsOptim = argsConfigOptim$numCallsOptim,
                control = argsConfigOptim$control,
                verbose = debug)

            } else {
              matParInitRunif <- PCMParamRandomVecParams(
                o = modelForFit,
                k = PCMNumTraits(modelForFit),
                R = PCMNumRegimes(modelForFit),
                n = argsConfigOptim$numRunifInitVecParams,
                argsPCMParamLowerLimit = argsConfigOptim$argsPCMParamLowerLimit,
                argsPCMParamUpperLimit = argsConfigOptim$argsPCMParamUpperLimit)
              matParInitGuess <- GuessInitVecParams(
                o = modelForFit,
                n = 1L,
                argsPCMParamLowerLimit = argsConfigOptim$argsPCMParamLowerLimit,
                argsPCMParamUpperLimit = argsConfigOptim$argsPCMParamUpperLimit,
                X = X, tree = tree, SE = SE,
                varyParams = FALSE)
              matParInitGuessVaryParams <- GuessInitVecParams(
                o = modelForFit,
                n = argsConfigOptim$numGuessInitVecParams,
                argsPCMParamLowerLimit = argsConfigOptim$argsPCMParamLowerLimit,
                argsPCMParamUpperLimit = argsConfigOptim$argsPCMParamUpperLimit,
                X = X, tree = tree, SE = SE,
                varyParams = TRUE)
              matParInitJitter <-
                jitterModelParams(
                  modelForFit,
                  argsPCMParamLowerLimit = argsConfigOptim$argsPCMParamLowerLimit,
                  argsPCMParamUpperLimit = argsConfigOptim$argsPCMParamUpperLimit,
                  numJitterRootRegimeFit = argsConfigOptim$numJitterRootRegimeFit,
                  sdJitterRootRegimeFit = argsConfigOptim$sdJitterRootRegimeFit,
                  numJitterAllRegimeFits = argsConfigOptim$numJitterAllRegimeFits,
                  sdJitterAllRegimeFits = argsConfigOptim$sdJitterAllRegimeFits,
                  verbose = debug)

              matParInit <- rbind(matParInitRunif,
                                  matParInitGuess,
                                  matParInitGuessVaryParams,
                                  matParInitJitter)

              if(nrow(fit) == 1L) {
                cat("Completing matParInit with vector from previous fit.\n")
                matParInit <- rbind(
                  matParInit, fit$fitVector[[1]][seq_len(ncol(matParInit))])
                print(matParInit[tail(nrow(matParInit)),])
              }

              fit <- PCMFit(
                X = XForFit, tree = treeForFit, model = modelForFit,
                SE = SEForFit,
                metaI = metaIFun, positiveValueGuard = positiveValueGuard,
                matParInit = matParInit,
                numCallsOptim = argsConfigOptim$numCallsOptim,
                control = argsConfigOptim$control,
                verbose = debug
              )
            }

            ll <- unname(logLik(fit))
            v_score = unname(scoreFun(fit))
            vec <- c(
              coef(fit),
              logLik = ll,
              df = attr(ll, "df"),
              nobs = attr(ll, "nobs"),
              score = v_score)

            dt.row <- data.table(
              hashCodeTree = hashCodes$hashCodeTree,
              hashCodeStartingNodesRegimesLabels =
                hashCodes$hashCodeStartingNodesRegimesLabels,
              hashCodeMapping = hashCodes$hashCodeMapping,
              treeEDExpression = EDExpression,
              startingNodesRegimesLabels = list(cladePartition),
              mapping = list(MatchModelMapping(modelMapping, modelTypes)),
              fitVector = list(unname(vec)),
              logLik = ll,
              df = attr(ll, "df"),
              nobs = attr(ll, "nobs"),
              score = v_score,
              duplicated = FALSE)
          }

          if(printFitVectorsToConsole) {
            cat(dt.row$treeEDExpression[[1]],
                ":",
                toString(unname(dt.row$fitVector[[1]])),
                "\n", sep="")
          }

          if(saveTempWorkerResults) {
            SaveTempWorkerResults(dt.row, prefixFiles)
          }
          dt.row
        }, silent=FALSE)

      } # end of nested foreach body

  CleanTemporaryFitFiles(filePrefix = prefixFiles)
  setkey(
    fits, hashCodeTree, hashCodeStartingNodesRegimesLabels, hashCodeMapping)
  fits
}
