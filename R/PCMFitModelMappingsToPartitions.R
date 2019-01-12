PCMFitModelMappingsToCladePartitions <- function(
  X, tree, modelTypes,
  SE = matrix(0.0, nrow(X), PCMTreeNumTips(tree)),

  listPartitions = NULL,
  listAllowedModelTypesIndices = NULL,

  scoreFun = AIC,

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

  fits <-
    foreach(
      cladePartition = listPartitions,
      EDExpression = EDExpressions,
      .combine = function(...) rbindlist(list(...)),
      .multicombine = TRUE,
      .inorder = FALSE,
      .packages = (.packages()) ) %:%

    foreach(
      modelMapping = PCMIteratorMapping2(
        mapping = unlist(
          sapply(
            listAllowedModelTypesIndices[as.character(cladePartition)],
            function(.) .[1])),
        modelTypes = seq_len(length(modelTypes)),
        allowedModelTypesIndices = lapply(
          listAllowedModelTypesIndices[as.character(cladePartition)],
          sort)),

      .combine=function(...) {
        CombineTaskResults(
          ...,
          filePrefix = prefixFiles,
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
            PCMTreeSetDefaultRegime(treeForFit, 1)
          } else {
            treeForFit <- tree
            PCMTreeSetRegimes(treeForFit, nodes = as.integer(cladePartition))
            XForFit <- X
            SEForFit <- SE
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
              # create an MixedGaussian model
              modelForFit <- do.call(MixedGaussian,
                                     c(list(k = nrow(XForFit),
                                            modelTypes = modelTypes,
                                            mapping = modelMapping),
                                       argsMixedGaussian))

              matParInitRunif <- PCMParamRandomVecParams(
                o = modelForFit,
                k = PCMNumTraits(modelForFit),
                R = PCMNumRegimes(modelForFit),
                n = argsConfigOptim$numRunifInitVecParams,
                argsPCMParamLowerLimit = argsConfigOptim$argsPCMParamLowerLimit,
                argsPCMParamUpperLimit = argsConfigOptim$argsPCMParamUpperLimit)
              matParInitGuess <- GuessInitVecParams(
                o = modelForFit,
                k = PCMNumTraits(modelForFit),
                R = PCMNumRegimes(modelForFit),
                n = argsConfigOptim$numGuessInitVecParams,
                argsPCMParamLowerLimit = argsConfigOptim$argsPCMParamLowerLimit,
                argsPCMParamUpperLimit = argsConfigOptim$argsPCMParamUpperLimit,
                X = X, tree = tree, SE = SE,
                varyTheta = FALSE)
              matParInitGuessVaryTheta <- GuessInitVecParams(
                o = modelForFit,
                k = PCMNumTraits(modelForFit),
                R = PCMNumRegimes(modelForFit),
                n = argsConfigOptim$numGuessInitVecParams,
                argsPCMParamLowerLimit = argsConfigOptim$argsPCMParamLowerLimit,
                argsPCMParamUpperLimit = argsConfigOptim$argsPCMParamUpperLimit,
                X = X, tree = tree, SE = SE,
                varyTheta = TRUE)

              matParInit <- rbind(matParInitRunif,
                                  matParInitGuess,
                                  matParInitGuessVaryTheta)

              fit <- PCMFit(
                X = XForFit, tree = treeForFit, model = modelForFit,
                SE = SEForFit,
                metaI = metaIFun, positiveValueGuard = positiveValueGuard,
                matParInit = matParInit,
                numCallsOptim = argsConfigOptim$numCallsOptim,
                control = argsConfigOptim$control,
                verbose = debug)

            } else {
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

              matParInitRunif <- PCMParamRandomVecParams(
                o = modelForFit,
                k = PCMNumTraits(modelForFit),
                R = PCMNumRegimes(modelForFit),
                n = argsConfigOptim$numRunifInitVecParams,
                argsPCMParamLowerLimit = argsConfigOptim$argsPCMParamLowerLimit,
                argsPCMParamUpperLimit = argsConfigOptim$argsPCMParamUpperLimit)
              matParInitGuess <- GuessInitVecParams(
                o = modelForFit,
                k = PCMNumTraits(modelForFit),
                R = PCMNumRegimes(modelForFit),
                n = argsConfigOptim$numGuessInitVecParams,
                argsPCMParamLowerLimit = argsConfigOptim$argsPCMParamLowerLimit,
                argsPCMParamUpperLimit = argsConfigOptim$argsPCMParamUpperLimit,
                X = X, tree = tree, SE = SE,
                varyTheta = FALSE)
              matParInitGuessVaryTheta <- GuessInitVecParams(
                o = modelForFit,
                k = PCMNumTraits(modelForFit),
                R = PCMNumRegimes(modelForFit),
                n = argsConfigOptim$numGuessInitVecParams,
                argsPCMParamLowerLimit = argsConfigOptim$argsPCMParamLowerLimit,
                argsPCMParamUpperLimit = argsConfigOptim$argsPCMParamUpperLimit,
                X = X, tree = tree, SE = SE,
                varyTheta = TRUE)
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
                                  matParInitGuessVaryTheta,
                                  matParInitJitter)

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
              treeEDExpression = EDExpression,
              hashCodeTree = hashCodes$hashCodeTree,
              hashCodeStartingNodesRegimesLabels =
                hashCodes$hashCodeStartingNodesRegimesLabels,
              hashCodeMapping = hashCodes$hashCodeMapping,
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
                dt.row$hashCodeTree[[1]],
                ":",
                dt.row$hashCodeStartingNodesRegimes[[1]],
                ":",
                dt.row$hashCodeMapping[[1]],
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
  fits
}
