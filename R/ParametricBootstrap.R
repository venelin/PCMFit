#' @export
CollectBootstrapResults <- function(
  resultDir,
  prefixFiles,
  ids,
  inferredTree,
  inferredBackboneTree,
  epochs,
  naturalNodeLabelsInOriginalTree = NULL,
  minLength = 0.2,
  verbose = FALSE) {

  rbindlist(
    lapply(
      ids,
      function(id) {
        if(verbose && id %% 1L == 0) {
          cat("id:", id, "\n")
        }

        prefixFilesId <- paste0(prefixFiles, id)

        resultFile <- paste0(
          resultDir, "/", prefixFilesId, "/FinalResult_", prefixFilesId, ".RData")

        CollectBootstrapResultInternal(
          id = id,
          resultFile = resultFile,
          inferredTreeHD = inferredTree,
          inferredBackboneTreeHD = inferredBackboneTree,
          epochsHD = epochs,
          newNodeLabelsInBsInferredTree = naturalNodeLabelsInOriginalTree,
          minLength = minLength,
          verbose = verbose
        )
      }))
}

#' @importFrom PCMBase TruePositiveRate FalsePositiveRate PCMTreeDtNodes PCMTreeSetLabels PCMMean PCMVar PCMTreeNearestNodesToEpoch PCMTreeGetPartsForNodes PCMTreeGetLabels PCMNumTraits PCMTreeGetPartNames PCMTreeMatrixNodesInSamePart
CollectBootstrapResultInternal <- function(
  id,
  resultFile,
  inferredTreeHD,
  inferredBackboneTreeHD,
  epochsHD,
  newNodeLabelsInBsInferredTree = NULL,
  minLength = 0.2,
  verbose = FALSE) {

  # Prevent check problems by creating variables with no visible binding:
  median <- mapping <- regime <- endNode <- SEs <- fitMappings <- NULL

  if(file.exists(resultFile)) {
    # this should load an object called fitMappings
    load(resultFile)
    fit <- RetrieveBestFitScore(fitMappings)

    bsInferredModel <- fit$inferredModel
    bsInferredModelTypes <- attr(bsInferredModel, "modelTypes")
    bsInferredMapping <- attr(bsInferredModel, "mapping")

    bsInferredTreeHD <- attr(bsInferredModel, "tree")
    ll <- attr(bsInferredModel, "ll")
    sc <- attr(bsInferredModel, "score")

    if(!is.null(newNodeLabelsInBsInferredTree)) {
      PCMTreeSetLabels(bsInferredTreeHD, newNodeLabelsInBsInferredTree)
    }
    # set natural labels to the non-artificially inserted internal nodes

    clusterNodesLabels <- list(
      PCMTreeGetLabels(bsInferredTreeHD)[as.integer(fit$inferredRegimeNodes)])

    # add singleton nodes at the past epochsHD for which we measure the
    # statistics of interest.
    for(epoch in epochsHD) {
      bsInferredTreeHD <-
        PCMTreeInsertSingletonsAtEpoch(bsInferredTreeHD, epoch, minLength = minLength)
    }

    bsInferredTreeHDDtNodes <- PCMTreeDtNodes(bsInferredTreeHD)
    bsInferredTreeHDDtNodes[, mapping:=bsInferredMapping[regime]]
    setkey(bsInferredTreeHDDtNodes, endNode)


    bsInferredMeans <- try(PCMMean(
      bsInferredTreeHD, bsInferredModel, internal = TRUE), silent = TRUE)
    bsInferredVar <- try(PCMVar(
      bsInferredTreeHD, bsInferredModel, SE = SEs, internal = TRUE)$Wii,
      silent = TRUE)

    if(! (class(bsInferredMeans) == "try-error" ||
          class(bsInferredVar) == "try-error") ) {

      # Compute global cluster tpr and fpr for each epoch;
      # Compute statistics at the intersection points of each branch
      # in inferredBackboneTreeHD with each epoch.
      listStatsForEpochs <- lapply(epochsHD, function(epoch) {
        # should be the same nodes for inferredTreeHD and bsInferredTreeHD!
        nodesAtEpoch <- PCMTreeNearestNodesToEpoch(inferredTreeHD, epoch)
        partsNodesAtEpoch <- PCMTreeGetPartsForNodes(inferredTreeHD, nodesAtEpoch)
        bsNodesAtEpoch <- PCMTreeNearestNodesToEpoch(bsInferredTreeHD, epoch)

        if(!identical(nodesAtEpoch, bsNodesAtEpoch)) {
          stop("nodesAtEpoch differs from bsNodesAtEpoch!")
        }

        trueClustersAtEpoch <-
          PCMTreeMatrixNodesInSamePart(inferredTreeHD, nodesAtEpoch)
        predClustersAtEpoch <-
          PCMTreeMatrixNodesInSamePart(bsInferredTreeHD, nodesAtEpoch)

        statsForEpoch <- list(
          tprCluster = TruePositiveRate(
            predClustersAtEpoch, trueClustersAtEpoch),
          fprCluster = FalsePositiveRate(
            predClustersAtEpoch, trueClustersAtEpoch),
          R = length(
            unique(PCMTreeGetPartsForNodes(bsInferredTreeHD, nodesAtEpoch)))
        )

        backbNodesAtEpoch <- PCMTreeNearestNodesToEpoch(inferredBackboneTreeHD, epoch)
        partsBackbNodesAtEpoch <-
          PCMTreeGetPartsForNodes(inferredBackboneTreeHD, backbNodesAtEpoch)

        if(!setequal(unique(partsBackbNodesAtEpoch), unique(partsNodesAtEpoch)) ) {
          cat("partsBackbNodesAtEpoch should have the same parts as partsNodesAtEpoch:\n")
          cat("backbNodesAtEpoch: ")
          print(PCMTreeGetLabels(inferredBackboneTreeHD)[backbNodesAtEpoch])
          print(unique(partsBackbNodesAtEpoch))
          cat("nodesAtEpoch: ")
          print(PCMTreeGetLabels(inferredTreeHD)[nodesAtEpoch])
          print(unique(partsNodesAtEpoch))
          stop("stopping.")
        }

        statsForBackboneNodes = lapply(
          seq_along(backbNodesAtEpoch), function(bni) {
            part <- partsBackbNodesAtEpoch[bni]
            nodesInSamePart <- nodesAtEpoch[which(partsNodesAtEpoch == part)]

            predClustersForNodesInSamePartAtEpoch <-
              PCMTreeMatrixNodesInSamePart(bsInferredTreeHD, nodesInSamePart)
            tprClusterForNodesInSamePart = TruePositiveRate(
              predClustersForNodesInSamePartAtEpoch,
              rep(1, length(predClustersForNodesInSamePartAtEpoch)))

            # Among the nodesInSamePart calculate the frequencies of the model
            # types, according to bsInferredModel and bsInferredTreeHD
            bsModelTypeFreqs <- rep(0.0, length(bsInferredModelTypes))

            for(mt in bsInferredTreeHDDtNodes[list(nodesInSamePart), mapping])
              bsModelTypeFreqs[mt] <- bsModelTypeFreqs[mt] + 1

            bsModelTypeFreqs <- bsModelTypeFreqs / length(nodesInSamePart)

            if(verbose) {
              cat(
                "epoch: ", epoch,
                ", backBoneNode: ",
                PCMTreeGetLabels(inferredBackboneTreeHD)[backbNodesAtEpoch[bni]],
                ", part: ", part, ",\n",
                "nodesInSamePart: ",
                toString(PCMTreeGetLabels(inferredTreeHD)[nodesInSamePart]),
                "\n",
                "tprClusterForNodesInSamePart: ", tprClusterForNodesInSamePart,
                "\n",
                "bsModelTypeFreqs: ", toString(bsModelTypeFreqs), "\n")
            }

            traitMean <- apply(
              bsInferredMeans[, nodesInSamePart, drop = FALSE], 1, mean)
            traitMedi <- apply(
              bsInferredMeans[, nodesInSamePart, drop = FALSE], 1, median)

            k <- PCMNumTraits(bsInferredModel)

            slopes <- do.call(cbind, lapply(nodesInSamePart, function(i) {
              Sigma <- bsInferredVar[, (i-1)*k + seq_len(k), drop = FALSE]
              # Sigma[1,2]/Sigma[1,1]
              as.vector((Sigma / matrix(diag(Sigma), k, k))[upper.tri(Sigma)])
            }))

            intercepts <- do.call(cbind, lapply(nodesInSamePart, function(i) {
              # a k x k matrix of the means, e.g. for k = 3 this is:
              # mu1 mu2 mu3
              # mu1 mu2 mu3
              # mu1 mu2 mu3
              muMat <- matrix(bsInferredMeans[, i, drop = FALSE], k, k, byrow = TRUE)
              Sigma <- bsInferredVar[, (i-1)*k + seq_len(k), drop = FALSE]
              # Sigma[1,2]/Sigma[1,1]
              # a k x k matrix of regression slopes on its upper triangle,
              # e.g. for k = 3 this is:
              # 1 b12 b13
              #     1 b23
              #         1
              # Here b12 means the regression slope of x2 on x1.
              bMat <- Sigma / matrix(diag(Sigma), k, k)
              # a k x k matrix of intercepts on its upper triangle:
              # 0 a12 a13
              #     0 a23
              #         0
              aMat <- muMat - bMat * matrix(diag(muMat), k, k)
              as.vector(aMat[upper.tri(aMat)])
            }))

            if(verbose) {
              cat("intercepts: ", toString(round(intercepts, 2)), "\n")
              cat("slopes: ", toString(round(slopes, 2)), "\n")
            }
            slopeMean <- apply(slopes, 1, mean)
            slopeMedi <- apply(slopes, 1, median)

            interceptMean <- apply(intercepts, 1, mean)
            interceptMedi <- apply(intercepts, 1, median)

            list(
              tprClusterForNodesInSamePart = tprClusterForNodesInSamePart,
              modelTypeFreqs = bsModelTypeFreqs,
              traitMean = traitMean,
              traitMedian = traitMedi,
              slopeMean = slopeMean,
              slopeMedian = slopeMedi,
              interceptMean = interceptMean,
              interceptMedian = interceptMedi)

          })

        names(statsForBackboneNodes) <-
          PCMTreeGetLabels(inferredBackboneTreeHD)[backbNodesAtEpoch]
        statsForEpoch$statsForBackboneNodes = statsForBackboneNodes
        statsForEpoch
      })
      names(listStatsForEpochs) <- as.character(epochsHD)

      data.table(
        Id = id,
        clusterNodesLabels = clusterNodesLabels,
        mapping = list(fit$inferredMappingIdx),
        model = list(fit$inferredModel),
        logLik = ll,
        score = sc,
        listStatsForEpochs = list(listStatsForEpochs))
    } else {
      data.table(
        Id = id,
        clusterNodesLabels = list(NULL),
        mapping = list(NULL),
        model = list(NULL),
        logLik = NA_real_,
        score = NA_real_,
        listStatsForEpochs = list(NULL))
    }
  } else {
    data.table(
      Id = id,
      clusterNodesLabels = list(NULL),
      mapping = list(NULL),
      model = list(NULL),
      logLik = NA_real_,
      score = NA_real_,
      listStatsForEpochs = list(NULL))
  }
}

#' @export
ExtractTimeSeriesForTraitValue <- function(
  model,
  epochs = seq(0, max(PCMTreeNodeTimes(attr(model, "tree"))), length.out = 20),
  backboneTree = NULL,
  traitIndex = 1L,
  traitName = paste0("Trait_", traitIndex)) {

  # prevent no visible binding warnings during check:
  timeInterval <- timeNode <- partFinal <- NULL

  if(is.null(backboneTree)) {
    tree <- attr(model, "tree")
    if(is.null(tree)) {
      stop("ExtractTimeSeriesForTraitValue: Both, backboneTree and attr(model, 'tree') are NULL.")
    } else {
      len <- max(PCMTreeNodeTimes(tree, tipsOnly = TRUE))
      epochs <- epochs[epochs >= 0 & epochs <= len]

      minLen <- (max(epochs) - min(epochs)) / length(epochs) / 10
      for(epoch in epochs) {
        tree <- PCMTreeInsertSingletonsAtEpoch(tree, epoch, minLength = minLen)
      }

      backboneTree <- PCMTreeBackbonePartition(tree)
    }
  }

  mI <- PCMInfo(X = NULL, tree = backboneTree, model = model)
  inferredMeans <- PCMMean(
    backboneTree, model, metaI = mI,
    internal = TRUE)

  nodeTimesBackbone <-
    structure(PCMTreeNodeTimes(backboneTree),
              names = PCMTreeGetLabels(backboneTree))

  res <- rbindlist(lapply(
    PCMTreeGetPartNames(backboneTree), function(part) {
      treeBackbonePart <-
        PCMTreeBackbonePartition(backboneTree, partsToKeep = part)

      rbindlist(lapply(epochs, function(epoch) {
        nodeAtEpoch <- PCMTreeNearestNodesToEpoch(treeBackbonePart, epoch)
        nodeAtEpochLab <- PCMTreeGetLabels(treeBackbonePart)[nodeAtEpoch]

        i <-  match(nodeAtEpochLab, PCMTreeGetLabels(backboneTree))

        mu <- inferredMeans[traitIndex, i]

        data.table(
          partFinal = part,
          regimeFinal = PCMTreeGetPartRegimes(treeBackbonePart)[part],
          node = nodeAtEpoch, nodeLab = nodeAtEpochLab,
          partNode = PCMTreeGetPartsForNodes(treeBackbonePart, nodeAtEpoch),
          regimeNode = PCMTreeGetPartRegimes(treeBackbonePart)[
            PCMTreeGetPartsForNodes(treeBackbonePart, nodeAtEpoch)],
          timeNode = nodeTimesBackbone[nodeAtEpochLab],
          value = mu)
      }))
    }))
  columnNamesWithoutTraitName <- c("value")

  setnames(
    res,
    old = columnNamesWithoutTraitName,
    new = paste0(traitName, ".", columnNamesWithoutTraitName))

  res[
    , timeInterval:=c(0, timeNode[2:.N]-timeNode[1:(.N-1L)]),
    by = partFinal]

  res
}

#' @export
#' @importFrom PCMBase PCMTreeNodeTimes PCMTreeInsertSingletonsAtEpoch PCMTreeBackbonePartition PCMTreeGetPartRegimes
ExtractTimeSeriesForTraitRegression <- function(
  model,
  epochs = seq(0, max(PCMTreeNodeTimes(attr(model, "tree"))), length.out = 20),
  backboneTree = NULL,
  traitIndexX = 1L,
  traitIndexY = 2L,
  traitNameX = paste0("Trait_", traitIndexX),
  traitNameY = paste0("Trait_", traitIndexY) ) {

  # prevent no visible binding warnings during check:
  timeInterval <- timeNode <- partFinal <- NULL
  if(is.null(backboneTree)) {
    tree <- attr(model, "tree")
    if(is.null(tree)) {
      stop("ExtractTimeSeriesForTraitValue: Both, backboneTree and attr(model, 'tree') are NULL.")
    } else {
      len <- max(PCMTreeNodeTimes(tree, tipsOnly = TRUE))
      epochs <- epochs[epochs >= 0 & epochs <= len]

      minLen <- (max(epochs) - min(epochs)) / length(epochs) / 10
      for(epoch in epochs) {
        tree <- PCMTreeInsertSingletonsAtEpoch(tree, epoch, minLength = minLen)
      }

      backboneTree <- PCMTreeBackbonePartition(tree)
    }
  }

  nodeTimesBackbone <-
    structure(PCMTreeNodeTimes(backboneTree),
              names = PCMTreeGetLabels(backboneTree))

  k <- PCMNumTraits(model)
  mI <- PCMInfo(X = NULL, tree = backboneTree, model = model)

  inferredMeans <- PCMMean(backboneTree, model, metaI = mI, internal = TRUE)
  inferredVars <- PCMVar(backboneTree, model, metaI = mI, internal = TRUE)$Wii

  res <- rbindlist(lapply(
    PCMTreeGetPartNames(backboneTree), function(part) {
      treeBackbonePart <-
        PCMTreeBackbonePartition(backboneTree, partsToKeep = part)

      rbindlist(lapply(epochs, function(epoch) {
        nodeAtEpoch <- PCMTreeNearestNodesToEpoch(treeBackbonePart, epoch)
        nodeAtEpochLab <- PCMTreeGetLabels(treeBackbonePart)[nodeAtEpoch]

        i <-  match(nodeAtEpochLab, PCMTreeGetLabels(backboneTree))

        mu <- inferredMeans[, i]
        Sigma <- inferredVars[, (i-1)*k + (1:k)]
        slope <- Sigma[traitIndexX, traitIndexY] / Sigma[traitIndexX, traitIndexX]
        intercept <- mu[traitIndexY] - mu[traitIndexX]*slope

        data.table(
          partFinal = part,
          regimeFinal = PCMTreeGetPartRegimes(treeBackbonePart)[part],
          node = nodeAtEpoch, nodeLab = nodeAtEpochLab,
          partNode = PCMTreeGetPartsForNodes(treeBackbonePart, nodeAtEpoch),
          regimeNode = PCMTreeGetPartRegimes(treeBackbonePart)[
            PCMTreeGetPartsForNodes(treeBackbonePart, nodeAtEpoch)],
          timeNode = nodeTimesBackbone[nodeAtEpochLab],

          slope = slope,
          intercept = intercept )

      }))
    }))

  columnNamesWithoutTraitName <- c("slope", "intercept")

  setnames(
    res,
    old = columnNamesWithoutTraitName,
    new = paste0(traitNameY, ".on.", traitNameX, ".", columnNamesWithoutTraitName))

  res[
    , timeInterval:=c(0, timeNode[2:.N]-timeNode[1:(.N-1L)]),
    by = partFinal]

  res
}


#' @importFrom data.table setnames
#' @export
ExtractBSDataForTraitValue <- function(
  tableBSFits,
  inferredModel,
  epochs,
  inferredBackboneTree,
  traitIndex = 1L,
  traitName = paste0("Trait_", traitIndex)
) {

  # Prevent no visible binding warnings during check.
  listStatsForEpochs <- Id <- quantile <- timeInterval <- timeNode <-
    partFinal <- NULL

  mI <- PCMInfo(X = NULL, tree = inferredBackboneTree, model = inferredModel)
  inferredMeans <- PCMMean(
    inferredBackboneTree, inferredModel, metaI = mI, internal = TRUE)

  nodeTimesBackbone <-
    structure(PCMTreeNodeTimes(inferredBackboneTree),
              names = PCMTreeGetLabels(inferredBackboneTree))

  res <- rbindlist(lapply(
    PCMTreeGetPartNames(inferredBackboneTree), function(part) {
      treeBackbonePart <-
        PCMTreeBackbonePartition(inferredBackboneTree, partsToKeep = part)

      rbindlist(lapply(epochs, function(epoch) {
        nodeAtEpoch <- PCMTreeNearestNodesToEpoch(treeBackbonePart, epoch)
        nodeAtEpochLab <- PCMTreeGetLabels(treeBackbonePart)[nodeAtEpoch]

        i <-  match(nodeAtEpochLab, PCMTreeGetLabels(inferredBackboneTree))

        mu <- inferredMeans[traitIndex, i]

        bsValues <- tableBSFits[
          ,
          sapply(.I, function(i) {
            if( !is.null(listStatsForEpochs[[i]]) ) {
              listStatsForEpochs[[i]][[as.character(epoch)]]$
                statsForBackboneNodes[[nodeAtEpochLab]]$traitMean[traitIndex]
            } else {
              NA_real_
            }
          })]

        bsIds <- tableBSFits[, Id]

        quantileBsValues <- quantile(
          bsValues, probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1),  na.rm = TRUE)

        data.table(
          partFinal = part,
          regimeFinal = PCMTreeGetPartRegimes(treeBackbonePart)[part],
          node = nodeAtEpoch, nodeLab = nodeAtEpochLab,
          partNode = PCMTreeGetPartsForNodes(treeBackbonePart, nodeAtEpoch),
          regimeNode = PCMTreeGetPartRegimes(treeBackbonePart)[
            PCMTreeGetPartsForNodes(treeBackbonePart, nodeAtEpoch)],
          timeNode = nodeTimesBackbone[nodeAtEpochLab],

          value = mu,

          bs.0 = quantileBsValues["0%"],
          bs.025 = quantileBsValues["2.5%"],
          bs.25 = quantileBsValues["25%"],
          bs.5 = quantileBsValues["50%"],
          bs.75 = quantileBsValues["75%"],
          bs.975 = quantileBsValues["97.5%"],
          bs1 = quantileBsValues["100%"],
          bsValues = list(bsValues),

          bs.Ids = list(bsIds)  )


      }))
    }))
  columnNamesWithoutTraitName <-
    c("value",
      "bs.0",
      "bs.025",
      "bs.25",
      "bs.5",
      "bs.75",
      "bs.975",
      "bs1",
      "bsValues",
      "bs.Ids")

  setnames(
    res,
    old = columnNamesWithoutTraitName,
    new = paste0(traitName, ".", columnNamesWithoutTraitName))

  res[
    , timeInterval:=c(0, timeNode[2:.N]-timeNode[1:(.N-1L)]),
    by = partFinal]

  res
}

#' @export
ExtractBSDataForTraitRegression <- function(
  tableBSFits,
  inferredModel,
  epochs,
  inferredBackboneTree,
  traitIndexX = 1L,
  traitIndexY = 2L,
  traitNameX = paste0("Trait_", traitIndexX),
  traitNameY = paste0("Trait_", traitIndexY) ) {

  # prevent no visible binding warnings during check
  quantile <- Id <- timeInterval <- timeNode <- partFinal <-
    listStatsForEpochs <- NULL

  nodeTimesBackbone <-
    structure(PCMTreeNodeTimes(inferredBackboneTree),
              names = PCMTreeGetLabels(inferredBackboneTree))

  k <- PCMNumTraits(inferredModel)
  mI <- PCMInfo(X = NULL, tree = inferredBackboneTree, model = inferredModel)

  inferredMeans <- PCMMean(
    inferredBackboneTree, inferredModel, metaI = mI, internal = TRUE)
  inferredVars <- PCMVar(
    inferredBackboneTree, inferredModel, metaI = mI, internal = TRUE)$Wii

  res <- rbindlist(lapply(
    PCMTreeGetPartNames(inferredBackboneTree), function(part) {
      treeBackbonePart <-
        PCMTreeBackbonePartition(inferredBackboneTree, partsToKeep = part)

      rbindlist(lapply(epochs, function(epoch) {
        nodeAtEpoch <- PCMTreeNearestNodesToEpoch(treeBackbonePart, epoch)
        nodeAtEpochLab <- PCMTreeGetLabels(treeBackbonePart)[nodeAtEpoch]

        i <-  match(nodeAtEpochLab, PCMTreeGetLabels(inferredBackboneTree))

        mu <- inferredMeans[, i]
        Sigma <- inferredVars[, (i-1)*k + (1:k)]
        slope <- Sigma[traitIndexX, traitIndexY] / Sigma[traitIndexX, traitIndexX]
        intercept <- mu[traitIndexY] - mu[traitIndexX]*slope

        bsSlopes <- tableBSFits[
          ,
          sapply(.I, function(i) {
            if( !is.null(listStatsForEpochs[[i]]) ) {
              matrixSlopes <- matrix(NA_real_, k, k)

              bsSlopeVals <-
                listStatsForEpochs[[i]][[as.character(epoch)]]$
                statsForBackboneNodes[[nodeAtEpochLab]]$slopeMean

              if(sum(upper.tri(matrixSlopes)) == length(bsSlopeVals)) {
                matrixSlopes[upper.tri(matrixSlopes)] <- bsSlopeVals
              }

              matrixSlopes[traitIndexX, traitIndexY]
            } else {
              NA_real_
            }
          })]

        quantileBsSlopes <- quantile(
          bsSlopes, probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1),  na.rm = TRUE)

        bsIntercepts <- tableBSFits[
          ,
          sapply(.I, function(i) {
            if( !is.null(listStatsForEpochs[[i]]) ) {
              matrixIntercepts <- matrix(NA_real_, k, k)
              bsIntercVals <- listStatsForEpochs[[i]][[as.character(epoch)]]$
                statsForBackboneNodes[[nodeAtEpochLab]]$interceptMean
              if(sum(upper.tri(matrixIntercepts)) == length(bsIntercVals)) {
                matrixIntercepts[upper.tri(matrixIntercepts)] <- bsIntercVals
              }
              matrixIntercepts[traitIndexX, traitIndexY]
            } else {
              NA_real_
            }
          })]

        quantileBsIntercepts <- quantile(
          bsIntercepts, probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1),  na.rm = TRUE)

        bsIds <- tableBSFits[, Id]

        data.table(
          partFinal = part,
          regimeFinal = PCMTreeGetPartRegimes(treeBackbonePart)[part],
          node = nodeAtEpoch, nodeLab = nodeAtEpochLab,
          partNode = PCMTreeGetPartsForNodes(treeBackbonePart, nodeAtEpoch),
          regimeNode = PCMTreeGetPartRegimes(treeBackbonePart)[
            PCMTreeGetPartsForNodes(treeBackbonePart, nodeAtEpoch)],
          timeNode = nodeTimesBackbone[nodeAtEpochLab],

          slope = slope,
          intercept = intercept,

          bs.slope.0 = quantileBsSlopes["0%"],
          bs.slope.025 = quantileBsSlopes["2.5%"],
          bs.slope.25 = quantileBsSlopes["25%"],
          bs.slope.5 = quantileBsSlopes["50%"],
          bs.slope.75 = quantileBsSlopes["75%"],
          bs.slope.975 = quantileBsSlopes["97.5%"],
          bs.slope1 = quantileBsSlopes["100%"],
          bs.SlopeValues = list(bsSlopes),

          bs.intercept.0 = quantileBsIntercepts["0%"],
          bs.intercept.025 = quantileBsIntercepts["2.5%"],
          bs.intercept.25 = quantileBsIntercepts["25%"],
          bs.intercept.5 = quantileBsIntercepts["50%"],
          bs.intercept.75 = quantileBsIntercepts["75%"],
          bs.intercept.975 = quantileBsIntercepts["97.5%"],
          bs.intercept1 = quantileBsIntercepts["100%"],
          bs.InterceptValues = list(bsIntercepts),

          bs.Ids = list(bsIds) )

      }))
    }))

  columnNamesWithoutTraitName <-
    c("slope",
      "intercept",
      "bs.slope.0",
      "bs.slope.025",
      "bs.slope.25",
      "bs.slope.5",
      "bs.slope.75",
      "bs.slope.975",
      "bs.slope1",
      "bs.SlopeValues",
      "bs.intercept.0",
      "bs.intercept.025",
      "bs.intercept.25",
      "bs.intercept.5",
      "bs.intercept.75",
      "bs.intercept.975",
      "bs.intercept1",
      "bs.InterceptValues",
      "bs.Ids" )

  setnames(
    res,
    old = columnNamesWithoutTraitName,
    new = paste0(traitNameY, ".on.", traitNameX, ".", columnNamesWithoutTraitName))

  res[
    , timeInterval:=c(0, timeNode[2:.N]-timeNode[1:(.N-1L)]),
    by = partFinal]

  res
}

#' @export
ExtractBSDataModelTypeFreqs <- function(
  tableBSFits,
  inferredModel,
  inferredTree,
  epochs,
  inferredBackboneTree) {

  # prevent no visible binding warnings during check:
  mapping <- regime <- endNodeLab <- listStatsForEpochs <- Id <-
    timeInterval <- timeNode <- partFinal <- startTime <- endTime <- NULL

  nodeTimesBackbone <-
    structure(PCMTreeNodeTimes(inferredBackboneTree),
              names = PCMTreeGetLabels(inferredBackboneTree))

  modelTypes <- attr(inferredModel, "modelTypes")
  inferredMapping <- attr(inferredModel, "mapping")

  dtInferredNodes <- PCMTreeDtNodes(inferredTree)
  dtInferredNodes[, mapping:=inferredMapping[regime]]
  setkey(dtInferredNodes, endNodeLab)

  res <- rbindlist(lapply(
    PCMTreeGetPartNames(inferredBackboneTree), function(part) {
      treeBackbonePart <-
        PCMTreeBackbonePartition(inferredBackboneTree, partsToKeep = part)

      rbindlist(lapply(epochs, function(epoch) {

        nodeAtEpoch <- PCMTreeNearestNodesToEpoch(treeBackbonePart, epoch)
        nodeAtEpochLab <- PCMTreeGetLabels(treeBackbonePart)[nodeAtEpoch]

        i <-  match(nodeAtEpochLab, PCMTreeGetLabels(inferredBackboneTree))

        inferredModelType <- rep(0.0, length(modelTypes))
        inferredModelType[dtInferredNodes[list(nodeAtEpochLab), mapping]] <- 1.0

        bsModelTypeFreqs <- tableBSFits[
          ,
          do.call(cbind, lapply(.I, function(i) {
            if( !is.null(listStatsForEpochs[[i]]) ) {
              listStatsForEpochs[[i]][[as.character(epoch)]]$
                statsForBackboneNodes[[nodeAtEpochLab]]$modelTypeFreqs
            } else {
              rep(NA_real_, length(modelTypes))
            }
          }))]

        bsIds <- tableBSFits[, Id]

        data.table(
          partFinal = part,
          regimeFinal = PCMTreeGetPartRegimes(treeBackbonePart)[part],
          node = nodeAtEpoch, nodeLab = nodeAtEpochLab,
          partNode = PCMTreeGetPartsForNodes(treeBackbonePart, nodeAtEpoch),
          regimeNode = PCMTreeGetPartRegimes(treeBackbonePart)[
            PCMTreeGetPartsForNodes(treeBackbonePart, nodeAtEpoch)],
          timeNode = nodeTimesBackbone[nodeAtEpochLab],

          modelTypeIndex = seq_along(modelTypes),

          inferredModelType = inferredModelType,

          bsModelTypeFreqs = lapply(
            seq_along(modelTypes), function(modelTypeI) {
              bsModelTypeFreqs[modelTypeI,]
            }),
          bs.Ids = list(bsIds) )
      }))
    }))

  res[
    , timeInterval:=c(
      rep(0.0, length(modelTypes)),
      timeNode[(length(modelTypes) + 1L):.N] -
        timeNode[1:(.N-length(modelTypes))]),
    by = partFinal]

  res
}


#' @export
ExtractBSDataTprWithinParts <- function(
  tableBSFits,
  inferredModel,
  epochs,
  inferredBackboneTree) {

  # prevent no visible binding warnings during check:
  listStatsForEpochs <- Id <- timeInterval <- timeNode <- partFinal <- NULL

  nodeTimesBackbone <-
    structure(PCMTreeNodeTimes(inferredBackboneTree),
              names = PCMTreeGetLabels(inferredBackboneTree))

  res <- rbindlist(lapply(
    PCMTreeGetPartNames(inferredBackboneTree), function(part) {
      treeBackbonePart <-
        PCMTreeBackbonePartition(inferredBackboneTree, partsToKeep = part)

      rbindlist(lapply(epochs, function(epoch) {

        nodeAtEpoch <- PCMTreeNearestNodesToEpoch(treeBackbonePart, epoch)
        nodeAtEpochLab <- PCMTreeGetLabels(treeBackbonePart)[nodeAtEpoch]

        i <-  match(nodeAtEpochLab, PCMTreeGetLabels(inferredBackboneTree))

        bsTprWithinPart <- tableBSFits[
          ,
          sapply(.I, function(i) {
            if( !is.null(listStatsForEpochs[[i]]) ) {
              listStatsForEpochs[[i]][[as.character(epoch)]]$
                statsForBackboneNodes[[nodeAtEpochLab]]$
                tprClusterForNodesInSamePart
            } else {
              NA_real_
            }
          })]

        bsIds <- tableBSFits[, Id]

        data.table(
          partFinal = part,
          regimeFinal = PCMTreeGetPartRegimes(treeBackbonePart)[part],
          node = nodeAtEpoch, nodeLab = nodeAtEpochLab,
          partNode = PCMTreeGetPartsForNodes(treeBackbonePart, nodeAtEpoch),
          regimeNode = PCMTreeGetPartRegimes(treeBackbonePart)[
            PCMTreeGetPartsForNodes(treeBackbonePart, nodeAtEpoch)],
          timeNode = nodeTimesBackbone[nodeAtEpochLab],

          bsTprWithinPart = list(bsTprWithinPart),
          bs.Ids = list(bsIds) )
      }))
    }))

  res[
    , timeInterval:=c(0, timeNode[2:.N]-timeNode[1:(.N-1L)]),
    by = partFinal]

  res
}

#' @importFrom PCMBase PCMTreeDtNodes
#' @export
ExtractBSDataTprFprGlobal <- function(
  tableBSFits,
  inferredModel,
  inferredBackboneTree,
  epochs) {

  # prevent no visible binding warnings during check:
  listStatsForEpochs <- Id <- startTime <- endTime <- NULL

  dtNodesBackboneTree <- PCMTreeDtNodes(inferredBackboneTree)
  res <- rbindlist(lapply(epochs, function(epoch) {

    numParts = nrow(dtNodesBackboneTree[startTime < epoch & epoch <= endTime])

    bsNumParts <- tableBSFits[
      ,
      sapply(.I, function(i) {
        if( !is.null(listStatsForEpochs[[i]]) ) {
          listStatsForEpochs[[i]][[as.character(epoch)]]$R
        } else {
          NA_real_
        }
      })]

    bsTpr <- tableBSFits[
      ,
      sapply(.I, function(i) {
        if( !is.null(listStatsForEpochs[[i]]) ) {
          listStatsForEpochs[[i]][[as.character(epoch)]]$tprCluster
        } else {
          NA_real_
        }
      })]
    bsFpr <- tableBSFits[
      ,
      sapply(.I, function(i) {
        if( !is.null(listStatsForEpochs[[i]]) ) {
          listStatsForEpochs[[i]][[as.character(epoch)]]$fprCluster
        } else {
          NA_real_
        }
      })]

    bsIds <- tableBSFits[, Id]

    data.table(
      epoch = epoch,
      numParts = numParts,
      bsNumParts = bsNumParts,
      bsTpr = bsTpr,
      bsFpr = bsFpr,
      bsIndex = bsIds)
  }))

  res
}
