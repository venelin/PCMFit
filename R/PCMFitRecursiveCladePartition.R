PCMFitRecursiveCladePartition <- function(
  X, tree, modelTypes,
  SE = matrix(0.0, nrow(X), PCMTreeNumTips(tree),
              dimnames=list(NULL, as.character(1:PCMTreeNumTips(tree)))),

  listCladePartitions = NULL,
  listAllowedModelTypesIndices = NULL,

  scoreFun = AIC,

  generatePCMModelsFun = NULL,
  metaIFun = PCMInfo,

  positiveValueGuard = Inf,

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

  argsConfigOptimAndMCMC = NULL,

  numJitterRootRegimeFit = 100, sdJitterRootRegimeFit = 0.5,
  numJitterAllRegimeFits = 100, sdJitterAllRegimeFits = 0.5,

  preorderTree = NULL,
  tableAncestors = NULL,

  saveTempWorkerResults = TRUE,
  printFitVectorsToConsole = FALSE,

  doParallel = FALSE,

  verbose = TRUE,
  verbosePCMFit = FALSE,
  verboseComposeMixedGaussianFromFits = FALSE,
  verboseAdaptArgsConfigOptimAndMCMC = FALSE
) {

  # Copy all arguments into a list
  # We establish arguments$<argument-name> as a convention for accessing the
  # original argument value.
  arguments <- as.list(environment())

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

  nodeLabelsTree <- PCMTreeGetLabels(tree)

  treeRootInt <- PCMTreeNumTips(tree) + 1
  treeRootLabel <- nodeLabelsTree[treeRootInt]
  hashCodeEntireTree <- tableFits[
    sapply(startingNodesRegimesLabels, length) == 1 &
      sapply(startingNodesRegimesLabels,
             function(s) match(treeRootInt, s, nomatch = 0L)) == 1,
    hashCodeTree[1]]

  if(verbose) {
    cat("Step 2: Perform a recursive clade partition search of an optimal mixed model...\n")
  }

  # 2.1 identify the best clade-partition and mapping currently in tableFits.
  # This should be the best fit of one of the model types to the whole tree.
  queuePartitionRoots <- tableFits[list(hashCodeEntireTree)][
    sapply(startingNodesRegimesLabels, length) == 1L, {
      iMinScore <- which.min(score)
      list(
        level = 1L,
        node = treeRootInt,
        partitionParentNode = treeRootInt,
        hashCodeBestPartitionInitial = hashCodeStartingNodesRegimesLabels[iMinScore],
        hashCodeBestMappingInitial = hashCodeMapping[iMinScore],
        hashCodeBestPartitionLevel = hashCodeStartingNodesRegimesLabels[iMinScore],
        hashCodeBestMappingLevel = hashCodeMapping[iMinScore]
      )
    }]

  if(verbose) {
    cat("Step 2.1: Initiated queue of partition root nodes with the root of the original tree...\n")
  }

  # index of the head row in queuePartitionRoots
  headQPR <- 1L

  mainLoopHistory <- list()

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

    bestScore <- tableFits[
      queuePartitionRoots[headQPR,
                          list(hashCodeEntireTree,
                               hashCodeBestPartitionLevel,
                               hashCodeBestMappingLevel)],
      score[[1]]]

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
      headQPR_Score = bestScore)

    if(verbose) {
      cat("Step 2.2: headQPR=", headQPR, ": bestPartition/mapping/score: ",
          toString(bestPartition), " / ",
          toString(match(bestMapping, modelTypes)),
          "(", toString(bestMapping), ") / ", bestScore, " ; \n")
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

    if(is.null(listCladePartitions)) {
      # 2.4. listCladePartitionsSubtree: Create a list of all possible
      # clade-partitions of subtree into clades not smaller than
      # minCladeSizes[partitionRootLevel].
      if(verbose) {
        cat("Step 2.4: Generating list of clade partitions for the subtree rooted at '",
            partitionRootLabel, "':\n")
      }

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
             minCladeSizeLevel > dtTipsPerRegime[node%in%c(
               partitionRootLabel, labelsSubtree[partitionSubtree]), min(N)]) {
            NULL
          } else {
            nodeLabelsTree[partNodes2]
          }
        })

      listCladePartitions <- listCladePartitions[!sapply(listCladePartitions, is.null)]
    } else {
      # at this point the elements in listCladePartititions have to be integers
      cladeRootsSubtree <- labelsSubtree[as.integer(unlist(listCladePartitions))]
      # now we convert them to characters
      listCladePartitions <- lapply(listCladePartitions, as.character)
    }

    if(verbose) {
      cat("Step 2.5: List of clade partitions of the tree:\n")
      cat(do.call(
        paste,
            c(lapply(listCladePartitions, function(p) {
              toString(p)
            }),
            list(sep="\n"))))
      cat("\n")
    }
    mainLoopHistoryEntry[["listCladePartitions"]] <- listCladePartitions


    # 2.6. Prepare "seeds" for the modelMappingIterators for each clade-partition
    if(verbose) {
      cat("Step 2.6: Preparing allowed-model-indices for the modelMappingIterators for each clade-partition\n")
    }
    PCMTreeSetRegimes(tree, as.integer(bestPartition))

    if(is.list(arguments$listAllowedModelTypesIndices)) {

      listAllowedModelTypesIndices <-
        arguments$listAllowedModelTypesIndices[union(
          bestPartition, cladeRootsSubtree)]

    } else if(is.null(arguments$listAllowedModelTypesIndices) ||
              arguments$listAllowedModelTypesIndices == "all") {

      nodes <- union(bestPartition, cladeRootsSubtree)
      listAllowedModelTypesIndices <-
        replicate(length(nodes), seq_along(modelTypes), simplify = FALSE)
      names(listAllowedModelTypesIndices) <- as.character(nodes)

    } else if(arguments$listAllowedModelTypesIndices %in%
              c("best-clade", "best-clade-2")) {

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
            sapply(startingNodesRegimesLabels,
                   function(s) match(label, s, nomatch = 0L)) == 1,
            match(mapping[[which.min(score)]][1], modelTypes)]


          if(arguments$listAllowedModelTypesIndices == "best-clade-2") {

            iLabel <- as.integer(label)
            # we need the if(), because PCMTreeGetRegimeForNode returns an empty
            # vector for the root-node
            iRegime <- if(iLabel == treeRootInt) {
              1
            } else {
              PCMTreeGetRegimeForNode(tree, iLabel)
            }
            bestCladePartitionMapping <- match(bestMapping[iRegime], modelTypes)

          } else {
            bestCladePartitionMapping <- bestCladeRootMapping
          }

          listAllowedModelTypesIndices[[label]] <-
            unique(c(bestCladePartitionMapping, bestCladeRootMapping))

          # if a previously inferred tableFits was supplied, it might contain
          # clade-fits for model-types other than
          # the supplied modelTypes in this fit. These result in NA entries in
          # allowedModelTypes-vectors. We remove them here:
          listAllowedModelTypesIndices[[label]] <-
            listAllowedModelTypesIndices[[label]][
              !is.na(listAllowedModelTypesIndices[[label]])]
        }
      }
    }
    if(verbose) {
      cat("listAllowedModelTypesIndices:\n",
          do.call(
            paste,
            c(lapply(seq_along(listAllowedModelTypesIndices), function(i) {
              paste0("`", names(listAllowedModelTypesIndices)[i], "`: ",
                     toString(listAllowedModelTypesIndices[[i]]))
            }),
            list(sep="\n"))))
      cat("\n")
    }

    mainLoopHistoryEntry[[
      "listAllowedModelTypesIndices"]] <- listAllowedModelTypesIndices

    # 2.7 Nested foreach over the clade-partitions and the allowed model-type
    # mappings
    argumentsFitsToTree <- arguments[
      intersect(
        names(arguments),
        names(as.list(args(PCMFitModelMappingsToCladePartitions))))]

    argumentsFitsToTree$X <- X
    argumentsFitsToTree$tree <- tree
    argumentsFitsToTree$modelTypes <- modelTypes
    argumentsFitsToTree$SE <- SE
    argumentsFitsToTree$listCladePartitions = listCladePartitions
    argumentsFitsToTree$listAllowedModelTypesIndices =
      listAllowedModelTypesIndices
    argumentsFitsToTree$fitMappingsPrev <- NULL
    argumentsFitsToTree$tableFitsPrev <- tableFits
    argumentsFitsToTree$modelTypesInTableFitsPrev <- modelTypesInTableFitsPrev

    argumentsFitsToTree$prefixFiles <- paste0(prefixFiles,
                                              partitionRootLevel,
                                              "_cladeparts_")

    argumentsFitsToTree$preorderTree <- preorderTree
    argumentsFitsToTree$tableAncestors <- tableAncestors

    fitsToTree <- do.call(PCMFitModelMappingsToCladePartitions,
                          argumentsFitsToTree)

    # update tableFits
    tableFits <- UpdateTableFits(tableFits, fitsToTree)

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
      if(fitsToTree[, min(score)] < bestScore &&
         partitionRootLevel < maxCladePartitionLevel ) {
        # a partitioning with a better score was found

        # put new partitioning nodes and the partition root-node in the queue
        # with an augmented level
        queuePartitionRootsNew <- fitsToTree[, {
          iMinScore <- which.min(score)
          list(
            level = partitionRootLevel + 1L,
            node = as.integer(c(partitionRootLabel,
                                setdiff(startingNodesRegimesLabels[[iMinScore]],
                                        bestPartition))),
            partitionParentNode = as.integer(partitionRootLabel),
            hashCodeBestPartitionInitial =
              hashCodeStartingNodesRegimesLabels[iMinScore],

            hashCodeBestMappingInitial = hashCodeMapping[iMinScore],
            hashCodeBestPartitionLevel = hashCodeStartingNodesRegimesLabels[iMinScore],
            hashCodeBestMappingLevel = hashCodeMapping[iMinScore]
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

        # append the new nodes to the queue and reset the key
        queuePartitionRoots <- rbindlist(list(queuePartitionRoots,
                                              queuePartitionRootsNew))
      }
    }
  } # end of main loop: while(headQPR <= nrow(queuePartitionRoots))
  res <- list(
    hashCodeEntireTree = hashCodeEntireTree,
    tableFits = tableFits,
    queuePartitionRoots = queuePartitionRoots,
    mainLoopHistory = mainLoopHistory)

}
