#' Fit regime-assignments to (sub-)trees in a tree with different assigned model
#' types to each regime.
#'
#' @description This function performs multiple model fits of mixed regime models
#' (MixedGaussian) mapping different model-types (e.g. BM and OU) to different
#' regimes in a tree and testing different regime assignments to the branches in
#' the tree.
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

  colnames(X) <- colnames(SE) <- as.character(1:PCMTreeNumTips(tree))

  tableFits <- InitTableFits(modelTypes,
                             fitMappingsPrev,
                             tableFitsPrev,
                             modelTypesInTableFitsPrev,
                             verbose = verbose)

  modelTypesInTableFits <- attr(tableFits, "modelTypes")

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
    fitsToClades <- PCMFitModelMappingsToClades(
      X = X, tree = tree, modelTypes = modelTypes,
      SE = SE,

      generatePCMModelsFun = generatePCMModelsFun,
      metaIFun = metaIFun, positiveValueGuard = positiveValueGuard,

      lik = lik, prior = prior, input.data = input.data, config = config,

      fitMappingsPrev = NULL,
      tableFitsPrev = tableFits,
      modelTypesInTableFitsPrev = modelTypesInTableFitsPrev,

      skipFitWhenFoundInTableFits = skipFitWhenFoundInTableFits,

      prefixFiles = prefixFiles,

      minCladeSize = MIN_CLADE_SIZE,

      listPCMOptions = listPCMOptions,

      argsMixedGaussian = argsMixedGaussian,
      argsPCMParamLowerLimit = argsPCMParamLowerLimit,
      argsPCMParamUpperLimit = argsPCMParamUpperLimit,
      argsPCMParamLoadOrStore = argsPCMParamLoadOrStore,

      argsConfigOptimAndMCMC = argsConfigOptimAndMCMC1,

      preorderTree = preorderTree,
      tableAncestors = tableAncestors,

      printFitVectorsToConsole = printFitVectorsToConsole,

      doParallel = doParallel,

      verbose = verbose,
      verbosePCMFit = verbosePCMFit
    )

    # update tableFits with the entries in fitsToClades which were not already there.
    if(nrow(tableFits)>0) {
      tableFits <- rbindlist(list(tableFits, fitsToClades))
    } else {
      tableFits <- fitsToClades
    }
    tableFits <- tableFits[duplicated == FALSE]
    setkey(tableFits, hashCodeTree,hashCodeStartingNodesRegimesLabels,hashCodeMapping)

    SaveCurrentResults(list(tableFits = tableFits),
                       filePrefix = prefixFiles)

    # 2. Perform fits to clade-partitions with different model mappings
    # we need these variables throughout this step

    nodeLabelsTree <- PCMTreeGetLabels(tree)

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

    # Main loop:
    while(headQPR <= nrow(queuePartitionRoots)) {
      # 2.2. pop the first partition root from the queue
      partitionRootNode <- queuePartitionRoots[headQPR, node]
      partitionRootLabel <- nodeLabelsTree[partitionRootNode]
      partitionRootLevel <- queuePartitionRoots[headQPR, level]

      hashCodeBestPartition <-
        queuePartitionRoots[headQPR, hashCodeBestPartitionLevel]

      bestPartition <- tableFits[
        queuePartitionRoots[headQPR,
                            list(hashCodeEntireTree,
                                 hashCodeBestPartitionLevel,
                                 hashCodeBestMappingLevel)],
        startingNodesRegimesLabels[[1]]]

      bestMapping <- tableFits[
        queuePartitionRoots[headQPR,
                            list(hashCodeEntireTree,
                                 hashCodeBestPartitionLevel,
                                 hashCodeBestMappingLevel)],
        mapping[[1]]]

      bestAIC <- tableFits[
        queuePartitionRoots[headQPR,
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
          queuePartitionRoots[-(1:headQPR),
                              list(level, node, partitionParentNode)],
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

        cat("Performing clade-partitioning at partitionRootLabel=",
            partitionRootLabel, "; partitionRootLevel=", partitionRootLevel,
            "\n")
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

          # if a previously inferred tableFits was supplied, it might contain
          # clade-fits for model-types other than
          # the supplied modelTypes in this fit. These result in NA entries in
          # allowedModelTypes-vectors. We remove them here:
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

      # 2.7 Nested foreach over the clade-partitions and the allowed model-type
      # mappings
      fitsToTree <- PCMFitModelMappingsToCladePartitions(
        X = X, tree = tree, modelTypes = modelTypes,
        SE = SE,

        listCladePartitions = listCladePartitions,
        listAllowedModelTypesIndices = listAllowedModelTypesIndices,

        generatePCMModelsFun = generatePCMModelsFun,
        metaIFun = metaIFun, positiveValueGuard = positiveValueGuard,

        lik = lik, prior = prior, input.data = input.data, config = config,

        fitMappingsPrev = NULL,
        tableFitsPrev = tableFits,
        modelTypesInTableFitsPrev = modelTypesInTableFitsPrev,

        skipFitWhenFoundInTableFits = skipFitWhenFoundInTableFits,

        prefixFiles = paste0(prefixFiles, partitionRootLevel, "_cladeparts_"),

        listPCMOptions = listPCMOptions,

        argsMixedGaussian = argsMixedGaussian,
        argsPCMParamLowerLimit = argsPCMParamLowerLimit,
        argsPCMParamUpperLimit = argsPCMParamUpperLimit,
        argsPCMParamLoadOrStore = argsPCMParamLoadOrStore,

        argsConfigOptimAndMCMC = argsConfigOptimAndMCMC2,

        numJitterRootRegimeFit = numJitterRootRegimeFit,
        sdJitterRootRegimeFit = sdJitterRootRegimeFit,
        numJitterAllRegimeFits = numJitterAllRegimeFits,
        sdJitterAllRegimeFits = sdJitterAllRegimeFits,

        preorderTree = preorderTree,
        tableAncestors = tableAncestors,

        printFitVectorsToConsole = printFitVectorsToConsole,

        doParallel = doParallel,

        verbose = verbose,
        verbosePCMFit = verbosePCMFit,
        verboseComposeMixedGaussianFromFits = verboseComposeMixedGaussianFromFits,
        verboseAdaptArgsConfigOptimAndMCMC = verboseAdaptArgsConfigOptimAndMCMC
      )

      if(nrow(tableFits) > 0) {
        tableFits <- rbindlist(list(tableFits, fitsToTree))
      } else {
        tableFits <- fitsToTree
      }
      tableFits <- tableFits[duplicated == FALSE]
      setkey(tableFits,
             hashCodeTree,hashCodeStartingNodesRegimesLabels,hashCodeMapping)

      # add new entry to the main loop history
      mainLoopHistory[[length(mainLoopHistory) + 1]] <- mainLoopHistoryEntry

      # Save the current results to a file "CurrentResults<filePrefix>.RData"
      SaveCurrentResults(list(tableFits = tableFits,
                              mainLoopHistory = mainLoopHistory,
                              queuePartitionRoots = queuePartitionRoots),
                         filePrefix = prefixFiles)

      # 2.8 Identify the best partition and mapping after the partitioning using
      # fitsToTree

      if( length(listCladePartitions) > 0 ) {
        if(fitsToTree[, min(aic)] < bestAIC &&
           partitionRootLevel < maxCladePartitionLevel ) {
          # a partitioning with a better AIC was found

          # put new partitioning nodes and the partition root-node in the queue
          # with an augmented level
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

          # update the best partition and mapping for the current and all
          # upcoming entries in the queue having
          # hashCodeBestPartitionLevel == hashCodeBestPartition with the current
          # level (partitionRootLevel) doing it in this way should facilitate
          # future inclusion of sub-optimal partitions in the queue
          hashCodeBestPartitionLevelNew <-
            queuePartitionRootsNew[, hashCodeBestPartitionInitial[1]]
          hashCodeBestMappingLevelNew <-
            queuePartitionRootsNew[, hashCodeBestMappingInitial[1]]
          queuePartitionRoots[
            c(rep(FALSE, headQPR - 1), rep(TRUE, .N - headQPR + 1)) &
              hashCodeBestPartitionLevel == hashCodeBestPartition,
            (c("hashCodeBestPartitionLevel", "hashCodeBestMappingLevel")) :=
              list(hashCodeBestPartitionLevelNew, hashCodeBestMappingLevelNew)]

          # previously the update was done only for the nodes at the current
          # level:
          # queuePartitionRoots[
          #   list(partitionRootLevel),
          #   hashCodeBestPartitionLevel:=hashCodeBestPartitionLevelNew]
          # queuePartitionRoots[
          #      list(partitionRootLevel),
          #      hashCodeBestMappingLevel:=hashCodeBestMappingLevelNew]

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
