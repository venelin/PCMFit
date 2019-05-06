#' Fit regime-assignments to (sub-)trees in a tree with different assigned model
#' types to each regime.
#'
#' @description This function performs multiple model fits of mixed regime models
#' (MixedGaussian) mapping different model-types (e.g. BM and OU) to different
#' regimes (colors) in a tree and testing different regime assignments to the
#' branches in the tree.
#' @importFrom foreach foreach when %do% %dopar% %:%
#' @importFrom data.table data.table rbindlist is.data.table setkey :=
#' @importFrom PCMBase PCMTree PCMTreeSetLabels PCMTreeSetPartition PCMTreeEvalNestedEDxOnTree PCMTreeNumTips PCMTreeListCladePartitions PCMTreeListAllPartitions PCMTreeToString MixedGaussian PCMOptions PCMTreeTableAncestors PCMTreeSplitAtNode PCMGetVecParamsRegimesAndModels MGPMDefaultModelTypes PCMGenerateModelTypes is.Transformable
#' @importFrom stats logLik coef AIC
#' @return an S3 object of class PCMFitModelMappings.
#'
#' @export
PCMFitMixed <- function(
  X, tree,

  modelTypes = MGPMDefaultModelTypes(),
  subModels = c(B = 'A', C = 'A', D = 'B', E = 'D', F = 'E'),

  argsMixedGaussian = Args_MixedGaussian_MGPMDefaultModelTypes(),

  SE = matrix(0.0, nrow(X), PCMTreeNumTips(tree)),

  generatePCMModelsFun = PCMGenerateModelTypes,

  metaIFun = PCMInfo, positiveValueGuard = Inf,

  scoreFun = AIC,

  fitMappingsPrev = NULL,
  tableFitsPrev = fitMappingsPrev$tableFits,
  modelTypesInTableFitsPrev = NULL,

  listPartitions = NULL,
  minCladeSizes = 20L,

  maxCladePartitionLevel = if(is.null(listPartitions)) 8L else 1L,
  maxNumNodesPerCladePartition = 1L,

  listAllowedModelTypesIndices = c("best-clade-2", "best-clade", "all"),

  argsConfigOptim1 = DefaultArgsConfigOptim(numCallsOptim = 10),
  argsConfigOptim2 = DefaultArgsConfigOptim(numCallsOptim = 4),
  argsConfigOptim3 = DefaultArgsConfigOptim(numCallsOptim = 10),

  maxNumRoundRobins = 0,
  maxNumPartitionsInRoundRobins = 2,

  listPCMOptions = PCMOptions(),

  skipFitWhenFoundInTableFits = TRUE,

  doParallel = FALSE,

  prefixFiles = "fits_",

  saveTempWorkerResults = TRUE,
  printFitVectorsToConsole = FALSE,

  verbose = TRUE,
  debug = FALSE
) {

  if( !is.null(listPartitions) ) {
    maxCladePartitionLevel = 1L
  }

  # Copy all arguments into a list
  # We establish arguments$<argument-name> as a convention for accessing the
  # original argument value.
  arguments <- as.list(environment())

  tree <- PCMTree(tree)
  PCMTreeSetLabels(tree)
  PCMTreeSetPartition(tree)

  colnames(X) <- colnames(SE) <- as.character(seq_len(PCMTreeNumTips(tree)))

  tableFits <- InitTableFits(modelTypes,
                             fitMappingsPrev,
                             tableFitsPrev,
                             modelTypesInTableFitsPrev,
                             verbose = verbose)

  modelTypesInTableFits <- attr(tableFits, "modelTypes")

  if(length(minCladeSizes) < maxCladePartitionLevel) {
    minCladeSizes <- c(
      rep(as.integer(NA), maxCladePartitionLevel - length(minCladeSizes)),
      minCladeSizes)
  }
  MIN_CLADE_SIZE <- min(minCladeSizes, na.rm = TRUE)

  if(PCMTreeNumTips(tree) > MIN_CLADE_SIZE) {

    preorderTree <- PCMTreePreorder(tree)
    tableAncestors <- PCMTreeTableAncestors(tree, preorder = preorderTree)

    # 1. (fitsToClades) Perform a fit of each model-type to each clade
    if(is.null(arguments$listPartitions) || arguments$listPartitions == "all") {
      cladeRoots <- c(PCMTreeNumTips(tree) + 1,
                      unlist(PCMTreeListCladePartitions(
                        tree = tree,
                        nNodes = 1,
                        minCladeSize = MIN_CLADE_SIZE,
                        tableAncestors = tableAncestors)))
    } else {
      cladeRoots = unique(c(PCMTreeNumTips(tree) + 1,
                            unlist(listPartitions)))
    }

    # prepare a list of allowed model type index vectors for the Fits To Clades
    listAllowedModelTypesIndicesFTC <-
      replicate(length(cladeRoots), seq_along(modelTypes), simplify = FALSE)
    names(listAllowedModelTypesIndicesFTC) <- as.character(cladeRoots)

    if(verbose) {
      cat("Step 1 (", Sys.time() ,"): Performing fits on", length(cladeRoots),
          " clades; ",
          sum(sapply(listAllowedModelTypesIndicesFTC[as.character(cladeRoots)],
                     length)), " model mappings altogether...\n",
          "Step 1.1 (", Sys.time() ,
          "): Fitting models independently from random starting locations...\n")
    }

    argumentsFitsToClades <-
      arguments[
        intersect(
          names(arguments),
          names(as.list(args(PCMFitModelMappingsToCladePartitions))))]

    argumentsFitsToClades$X <- X
    argumentsFitsToClades$tree <- tree
    argumentsFitsToClades$modelTypes <- modelTypes
    argumentsFitsToClades$SE <- SE
    argumentsFitsToClades$listPartitions <- as.list(cladeRoots)
    argumentsFitsToClades$listAllowedModelTypesIndices <-
      listAllowedModelTypesIndicesFTC
    argumentsFitsToClades$fitClades <- TRUE
    argumentsFitsToClades$fitMappingsPrev <- NULL
    argumentsFitsToClades$tableFitsPrev <- tableFits
    argumentsFitsToClades$modelTypesInTableFitsPrev <- modelTypesInTableFitsPrev
    argumentsFitsToClades$argsConfigOptim <- argsConfigOptim1
    argumentsFitsToClades$preorderTree <- preorderTree
    argumentsFitsToClades$tableAncestors <- tableAncestors
    argumentsFitsToClades$prefixFiles <- paste0(prefixFiles, "_clades_")

    fitsToClades <- do.call(
      PCMFitModelMappingsToCladePartitions, argumentsFitsToClades)

    # Fix suboptimal fits, in which a sub-model of the fitted model got a higher
    # likelihood value.
    checkForBetterSubmodels <- TRUE
    checkForBetterSubmodelsIteration <- 0L
    while(checkForBetterSubmodels) {
      checkForBetterSubmodelsIteration <- checkForBetterSubmodelsIteration + 1L
      if(verbose) {
        cat(
          "Step 1.2, Iteration ", checkForBetterSubmodelsIteration,
          "(", Sys.time() ,"):",
          "Learning from sub-models, where the found max log-likelihood of a",
          "super-model was lower than the one of its sub-model...\n")
      }

      betterSubmodelFits <- LearnCladeFitsFromSubmodels(
        cladeFits = fitsToClades,
        modelTypes = modelTypes,
        subModels = subModels,
        argsMixedGaussian = argsMixedGaussian,
        metaIFun = metaIFun,
        scoreFun = scoreFun,
        X = X, tree = tree, SE = SE,
        verbose = verbose)

      if(nrow(betterSubmodelFits$cladeFitsNew) > 0L) {

        fitsToClades <- UpdateTableFits(
          fitsToClades, betterSubmodelFits$cladeFitsNew)

        argumentsFitsToClades$listPartitions <-
          betterSubmodelFits$listPartitions
        argumentsFitsToClades$listAllowedModelTypesIndices <-
          betterSubmodelFits$listAllowedModelTypesIndices
        argumentsFitsToClades$skipFitWhenFoundInTableFits <- FALSE
        argumentsFitsToClades$argsConfigOptim <-
          DefaultArgsConfigOptim(
            numRunifInitVecParams = 2L,
            numGuessInitVecParams = 2L,
            numJitterRootRegimeFit = 2L,
            numJitterAllRegimeFits = 2L,
            numCallsOptim = 1L)
        argumentsFitsToClades$tableFitsPrev <-
          betterSubmodelFits$cladeFitsNew

        fitsToCladesRerun <- do.call(
          PCMFitModelMappingsToCladePartitions, argumentsFitsToClades)

        fitsToClades <- UpdateTableFits(fitsToClades, fitsToCladesRerun)
      } else {
        checkForBetterSubmodels <- FALSE
      }
    }


    # update tableFits with the entries in fitsToClades
    tableFits <- UpdateTableFits(tableFits, fitsToClades)

    SaveCurrentResults(list(tableFits = fitsToClades), filePrefix = prefixFiles)

    # 2. Perform fits to clade-partitions with different model mappings
    # we need these variables throughout this step
    if(!is.list(arguments$listAllowedModelTypesIndices)) {
      # by default listAllowedModelTypesIndices is a character vector listing
      # all possible values. Here, we retain only the first value for the
      # subsequent call to PCMFitRecursiveCladePartition.
      arguments$listAllowedModelTypesIndices <-
        arguments$listAllowedModelTypesIndices[1]
    }
    argumentsStep2 <- arguments[intersect(
      names(arguments), names(as.list(args(PCMFitRecursiveCladePartition))))]

    argumentsStep2$X <- X
    argumentsStep2$tree <- tree
    argumentsStep2$modelTypes <- modelTypes
    argumentsStep2$SE <- SE
    argumentsStep2$fitMappingsPrev <- NULL
    argumentsStep2$tableFitsPrev <- tableFits
    argumentsStep2$modelTypesInTableFitsPrev <- modelTypesInTableFitsPrev
    argumentsStep2$argsConfigOptim <- argsConfigOptim2
    argumentsStep2$preorderTree <- preorderTree
    argumentsStep2$tableAncestors <- tableAncestors

    resultStep2 <- do.call(PCMFitRecursiveCladePartition, argumentsStep2)

    fitsToTree <- rbindlist(list(
      fitsToClades[list(hashCodeTree = resultStep2$hashCodeEntireTree)],
      resultStep2$fitsToTree))

    # update tableFits with the entries in fitsToTree
    tableFits <- UpdateTableFits(tableFits, fitsToTree)

    SaveCurrentResults(list(tableFits = tableFits), filePrefix = prefixFiles)

    # A table with the best fit for each of the top
    # maxNumPartitionsInRoundRobins partitions in resultStep2$fitsToTree
    tableFitsRRInit <- fitsToTree[
      ,
      .SD[which.min(score)],
      keyby = hashCodeStartingNodesRegimesLabels][
        order(score)][
          seq_len(min(maxNumPartitionsInRoundRobins, .N)), .SD,
          keyby = hashCodeStartingNodesRegimesLabels]

    tableFitsRR <- NULL

    # Step 3. Round robin : This is an optional step controlled by the argument
    # maxNumRoundRobins, which is 0 by default.
    if(maxNumRoundRobins > 0) {

      tableFitsRR <- copy(tableFitsRRInit)

      canImprove <- rep(TRUE, nrow(tableFitsRR))

      if(verbose) {
        cat("Step 3 (", Sys.time() ,"): Performing up to", maxNumRoundRobins,
            "round robin iterations; initial selected partitions/mappings:\n")
        print(tableFitsRR[, list(
          hashCodeStartingNodesRegimesLabels,
          startingNodesRegimesLabels,
          mapping = lapply(mapping, function(m) {
            names(modelTypes)[match(m, modelTypes)]
          }),
          logLik,
          R = sapply(startingNodesRegimesLabels, length),
          df, score,
          canImprove = canImprove)])
      }

      iRR <- 1L
      while(iRR < maxNumRoundRobins) {
        dtOldScore <- tableFitsRR[
          , list(oldScore = score), keyby = hashCodeStartingNodesRegimesLabels]
        partitionLengths <- tableFitsRR[
          , sapply(startingNodesRegimesLabels, length)]

        for( pos in seq_len(max(partitionLengths)) ) {

          # logical vector indicating the partitions in `partitions` not shorter
          # than this pos
          haveThisPos <- (pos <= partitionLengths)

          if(verbose) {
            cat(
              "> Step 3, iteration ", iRR, " (", Sys.time(), ")",
              "round robin loop for node position", pos, "; ",
              sum(canImprove & haveThisPos),
              "of the selected top partitions/mappings have this position and",
              "might improve their score. \n")
          }

          if(sum(canImprove & haveThisPos) > 0) {
            tableFitsRRForPos <- RetrieveFittedModelsFromFitVectors(
              fitMappings = NULL,
              tableFits = tableFitsRR[canImprove & haveThisPos],
              modelTypes = modelTypes,
              modelTypesNew = NULL,
              argsMixedGaussian = argsMixedGaussian,
              X = X,
              tree = tree,
              SE = SE,
              setAttributes = FALSE
            )

            partitions <- tableFitsRRForPos$startingNodesRegimesLabels
            mappings <- tableFitsRRForPos$mapping

            listHintModels <- tableFitsRRForPos$fittedModel

            listNamesInHintModels <- lapply(listHintModels, function(hm) {
              # all member names except "pos"
              ipos <- match(as.character(pos), names(hm))
              names(hm)[-ipos]
            })

            # create listAllowedModelTypesIndices for each partition in
            # partitions
            listAllowedModelTypesIndices <- lapply(
              seq_along(partitions),

              function(iPartition) {
                m <- lapply(mappings[[iPartition]], match, modelTypes)
                m[[pos]] <- seq_along(modelTypes)
                names(m) <- as.character(
                  partitions[[iPartition]])
                m
              })

            # Call PCMFitModelMappingsToPartitions
            argumentsRR <-
              arguments[
                intersect(
                  names(arguments),
                  names(as.list(args(PCMFitModelMappingsToCladePartitions))))]

            argumentsRR$X <- X
            argumentsRR$tree <- tree
            argumentsRR$modelTypes <- modelTypes
            argumentsRR$SE <- SE
            argumentsRR$listPartitions <- partitions
            argumentsRR$listAllowedModelTypesIndices <-
              listAllowedModelTypesIndices
            argumentsRR$fitClades <- FALSE
            argumentsRR$fitMappingsPrev <- NULL
            argumentsRR$tableFitsPrev <- fitsToClades
            argumentsRR$modelTypesInTableFitsPrev <- modelTypes
            argumentsRR$argsConfigOptim <- argsConfigOptim3
            argumentsRR$preorderTree <- preorderTree
            argumentsRR$tableAncestors <- tableAncestors
            argumentsRR$prefixFiles <- paste0(prefixFiles, "_rr_")
            argumentsRR$listHintModels <- listHintModels
            argumentsRR$listNamesInHintModels <- listNamesInHintModels


            newFitsForThisPos <- do.call(
              PCMFitModelMappingsToCladePartitions, argumentsRR)

            if(!is.data.table(newFitsForThisPos)) {
              cat("newFitsFortThisPos not a data.table, but is:\n")
              print(newFitsForThisPos)

              errorList <- list(argumentsRR = argumentsRR,
                                newFitsForThisPos = newFitsForThisPos)
              save(errorList, file="ListErrorObjects.RData")
            }

            # update fitsToTree
            fitsToTree <- UpdateTableFits(fitsToTree, newFitsForThisPos)

            # update tableFits with the entries in fitsToTree
            tableFits <- UpdateTableFits(tableFits, fitsToTree)

            SaveCurrentResults(
              list(tableFits = tableFits), filePrefix = prefixFiles)

            # update tableFitsRR with the new best fits
            # Here we do not use UpdateTableFits, because it uses another key
            tableFitsRR <- rbindlist(
              list(tableFitsRR, newFitsForThisPos),
              use.names = TRUE)
            # Keep the best mapping for each partition:
            tableFitsRR <- tableFitsRR[
              ,
              .SD[which.min(score)],
              keyby = hashCodeStartingNodesRegimesLabels]
          }

          if(verbose) {
            cat(
              "> Step 3, iteration ", iRR, " (", Sys.time(), ")",
              "round robin loop for node position", pos,
              ", top candidate partitions/mappings at the end of this loop:\n")
            print(tableFitsRR[, list(
              hashCodeStartingNodesRegimesLabels,
              startingNodesRegimesLabels,
              mapping = lapply(mapping, function(m) {
                names(modelTypes)[match(m, modelTypes)]
              }),
              logLik,
              R = sapply(startingNodesRegimesLabels, length),
              df, score,
              canImprove = canImprove)])
          }
        }

        canImprove <- tableFitsRR[dtOldScore, score < oldScore]

        if(sum(canImprove) == 0) {
          if(verbose) {
            cat("No scores improved during the last round robin iteration")
            if(iRR < maxNumRoundRobins) {
              cat("Exiting round robin before reaching", maxNumRoundRobins, "iterations.")
            }
          }
          break
        }

        iRR <- iRR + 1L
      }
    }
  } else {
    stop("ERR:04131:PCMFit:PCMFitMixed:PCMFitMixed:: the tree has fewer tips than the min clade-size in minCladeSizes (",
         MIN_CLADE_SIZE,
         "). Try with smaller minCladeSizes.")
  }

  arguments$tableFitsPrev <- arguments$fitMappingsPrev <- paste0(
    "To save space and avoid redundancy, tableFitsPrev and fitMappingsPrev\n",
    "are not stored in arguments. Most fits in tableFitsPrev are found\n",
    "in the member tableFits in the PCMFitModelMappings object. However some\n",
    "of these fits might have been replaced by equivalent model fits with a\n",
    "higher score encountered during this call to PCMFitMixed search.")

  resFitMappings <- list(
    arguments = arguments,
    options = PCMOptions(),
    tree = tree,
    X = X,
    SE = SE,
    hashCodeTree = resultStep2$hashCodeEntireTree,
    tableFits = tableFits,
    tableFitsThisSearchOnly = rbindlist(list(fitsToTree, fitsToClades), use.names = TRUE),
    queuePartitionRoots = resultStep2$queuePartitionRoots,
    mainLoopHistory = resultStep2$mainLoopHistory,
    tableFitsRRInit = tableFitsRRInit,
    tableFitsRR = tableFitsRR
  )
  class(resFitMappings) <- "PCMFitModelMappings"

  resFitMappings
}
