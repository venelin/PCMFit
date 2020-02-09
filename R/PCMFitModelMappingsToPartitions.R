#' Internally used cache environment.
#' @export
.cacheEnv <- new.env()

#' @importFrom utils tail
PCMFitModelMappingsToCladePartitions <- function(
  X, tree, modelTypes,
  treeVCVMat = NULL,
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

  cladePartition <- iCladePartition <- hintModel <-
    namesInHintModel <- EDExpression <- tableAncestorsTree <- NULL

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
      modelMapping = PCMIteratorMapping(
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
          # potentially big tableAncestors object
          if(!exists("tableAncestorsTree", .cacheEnv)) {
            if(verbose) {
              cat("Creating tableAncestorsTree to be used during fits\n")
            }
            assign("tableAncestorsTree",
                   PCMTreeTableAncestors(tree, preorderTree), .cacheEnv)
          }

          tableAncestorsTree <- get("tableAncestorsTree", envir = .cacheEnv)

          if(!exists("generatedPCMModels", .cacheEnv) &&
             !is.null(generatePCMModelsFun)) {
            if(verbose) {
              cat("Calling generatePCMModelsFun()...\n")
            }
            # this should generate PCMParentClasses and PCMSpecify functions
            # for all models in modelTypes
            generatePCMModelsFun()
            assign("generatedPCMModels", TRUE, .cacheEnv)
          }

          # don't want the names in modelMapping, the positions correspond to
          # cladePartition
          modelMapping <- unname(modelMapping)

          # set PCMBase options from parent process: necessary if this is
          # executed by a worker process from a cluster.
          options(listPCMOptions)

          if(fitClades) {
            k <- nrow(X)
            XInd <- rbind(X, as.double(seq_len(ncol(X))))
            treeSplit <- PCMTreeSplitAtNode(
              tree, as.integer(cladePartition[1]), tableAncestorsTree, XInd)
            treeForFit <- treeSplit$clade
            XForFit <- treeSplit$Xclade[1:k, , drop = FALSE]
            SEForFit <- if(is.matrix(SE)) {
              # SE is k x N matrix
              SE[, treeSplit$Xclade[k+1L, ], drop = FALSE]
            } else {
              # SE is a k x k x N cube
              SE[, , as.integer(treeSplit$Xclade[k+1L, ]), drop = FALSE]
            }
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
                X = X, tree = tree, treeVCVMat = treeVCVMat, SE = SE,
                tableAnc = tableAncestorsTree,
                varyParams = FALSE)
              matParInitGuessVaryParams <- do.call(
                rbind, lapply(
                  seq_len(argsConfigOptim$numGuessInitVecParams / 100),
                  function(i) {
                    GuessInitVecParams(
                      o = modelForFit,
                      n = 100,
                      argsPCMParamLowerLimit = argsConfigOptim$argsPCMParamLowerLimit,
                      argsPCMParamUpperLimit = argsConfigOptim$argsPCMParamUpperLimit,
                      X = X, tree = tree, treeVCVMat = treeVCVMat, SE = SE,
                      tableAnc = tableAncestorsTree,
                      varyParams = TRUE)
                  }))


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
                X = X, tree = tree, treeVCVMat = treeVCVMat, SE = SE,
                tableAnc = tableAncestorsTree,
                varyParams = FALSE)
              matParInitGuessVaryParams <- do.call(
                rbind, lapply(
                  seq_len(argsConfigOptim$numGuessInitVecParams / 100),
                  function(i) {
                    GuessInitVecParams(
                      o = modelForFit,
                      n = 100,
                      argsPCMParamLowerLimit = argsConfigOptim$argsPCMParamLowerLimit,
                      argsPCMParamUpperLimit = argsConfigOptim$argsPCMParamUpperLimit,
                      X = X, tree = tree, treeVCVMat = treeVCVMat, SE = SE,
                      tableAnc = tableAncestorsTree,
                      varyParams = TRUE)
                  }))
              matParInitJitter <-
                JitterModelParams(
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

            # prevent 'no visible binding' notes
            hashCodeTree <- hashCodeStartingNodesRegimesLabels <-
              hashCodeMapping <- NULL

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
            # define this internal function here, since it is invisible outside
            # the foreach body unless exported.
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
            SaveTempWorkerResults(dt.row, prefixFiles)
          }
          dt.row
        }, silent=FALSE)

      } # end of nested foreach body

  CleanTemporaryFitFiles(filePrefix = prefixFiles)
  if(is.data.table(fits)) {
    setkey(
      fits, hashCodeTree, hashCodeStartingNodesRegimesLabels, hashCodeMapping)
  }

  fits
}
