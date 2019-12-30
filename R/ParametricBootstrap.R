#' Collect the results from a parametric bootstrap MGPM fit
#' @param resultDir,prefixFiles character strings denoting the directory
#' and the file prefix used to locate the result .RData files. The way how the
#' result file-names are formed is specified by the argument
#' \code{exprResultFile}.
#' @param exprResultFile a character string evaluating to an R expression.
#' Default setting:
#' \code{
#' 'paste0(prefixFiles, id,
#'         "/Result_", prefixFiles, id, ".RData")'
#' }
#' In the above R-expression, \code{id} is any integer element from \code{ids}.
#' @param ids an integer vector denoting the ids of the bootstrap iterations.
#' @param inferredTree a PCMTree or a phylo object with set partition as
#' inferred during the ML or recursive clade partition search. This is the tree
#' used to simulate the bootstrap datasets.
#' @param epochs a double vector denoting in increasing order the time-points
#' at which ancestral mean and regression slope parameters should be calculated.
#' The inferredTree is to be 'sliced' by inserting singleton nodes at all
#' branch crossings with each epoch. By default the epochs are set at equal
#' intervals of length 10*minLength, via the following code:
#' \code{
#' seq(minLength, max(PCMTreeNodeTimes(inferredTree)), by = 10*minLength)}.
#' See also \code{minLength}.
#' @param minLength a positive double denoting the minimal time distance between
#' a node in the tree and an epoch in \code{epochs}, at which an additional
#' singleton node can be inserted. Default: 0.2.
#' @param nameFitObject a character string (default: 'fitMappings') denoting the
#' name of the PCMFitModelMappings object inside each result .RData file (see
#' also \code{resultDir} and \code{prefixFiles}).
#' @param keepBootstrapData logical indicating if the returned object should
#' contain the bootstrap simulated trait data.
#' @param verbose logical indicating if some debug or information messages
#' should be printed during the collection of the results.
#' @return a named list of S3 class 'ParametricBootstrapFits'.
#' @export
CollectBootstrapResults <- function(
  resultDir,
  prefixFiles,
  ids,
  inferredTree,
  epochs = seq(
    minLength, max(PCMTreeNodeTimes(inferredTree)), by = 10*minLength),
  minLength = 0.2,
  exprResultFile = 'paste0(prefixFiles, id, "/Result_", prefixFiles, id, ".RData")',
  nameFitObject = "fitMappings",
  keepBootstrapData = FALSE,
  verbose = FALSE) {

  inferredTreeHD <- inferredTree

  for(epoch in epochs) {
    inferredTreeHD <- PCMTreeInsertSingletonsAtEpoch(
      inferredTreeHD, epoch, minLength = minLength)
  }

  inferredBackboneTreeHD <- PCMTreeBackbonePartition(inferredTreeHD)
  naturalNodeLabelsInOriginalTree = PCMTreeGetLabels(inferredTree)

  tableBSFits <- rbindlist(
    lapply(
      ids,
      function(id) {
        if(verbose && id %% 1L == 0) {
          cat("id:", id, "\n")
        }

        resultFile <- paste0(resultDir, "/", eval(parse(text = exprResultFile)))

        CollectBootstrapResultInternal(
          id = id,
          resultFile = resultFile,
          inferredTreeHD = inferredTreeHD,
          inferredBackboneTreeHD = inferredBackboneTreeHD,
          epochsHD = epochs,
          newNodeLabelsInBsInferredTree = naturalNodeLabelsInOriginalTree,
          minLength = minLength,
          nameFitObject = nameFitObject,
          keepBootstrapData = keepBootstrapData,
          verbose = verbose)
      }))
  res <- list(
    inferredTree = inferredTree,
    inferredTreeHD = inferredTreeHD,
    epochs = epochs,
    minLength = minLength,
    tableBSFits = tableBSFits)
  class(res) <- "ParametricBootstrapFits"
  res
}

#' @importFrom PCMBase TruePositiveRate FalsePositiveRate PCMTreeDtNodes
#' PCMTreeSetLabels PCMMean PCMVar PCMTreeNearestNodesToEpoch
#' PCMTreeGetPartsForNodes PCMTreeGetLabels PCMNumTraits PCMTreeGetPartNames
#' PCMTreeMatrixNodesInSamePart PCMTreeInsertSingletonsAtEpoch
#' PCMTreeBackbonePartition
CollectBootstrapResultInternal <- function(
  id,
  resultFile,
  inferredTreeHD,
  inferredBackboneTreeHD,
  epochsHD,
  newNodeLabelsInBsInferredTree = NULL,
  minLength = 0.2,
  nameFitObject = "fitMappings",
  keepBootstrapData = FALSE,
  verbose = FALSE) {

  # Prevent check problems by creating variables with no visible binding:
  median <- mapping <- regime <- endNode <- fitMappings <- NULL

  if(file.exists(resultFile)) {
    # this should load an object called fitMappings
    load(resultFile)
    fit <- RetrieveBestFitScore(get(nameFitObject))

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
    clusterNodesLabels <- list(PCMTreeGetLabels(bsInferredTreeHD)[
      as.integer(fit$inferredRegimeNodes)])

    # add singleton nodes at the past epochsHD for which we measure the
    # statistics of interest.
    for(epoch in epochsHD) {
      bsInferredTreeHD <- PCMTreeInsertSingletonsAtEpoch(
        bsInferredTreeHD, epoch, minLength = minLength)
    }

    bsInferredTreeHDDtNodes <- PCMTreeDtNodes(bsInferredTreeHD)
    bsInferredTreeHDDtNodes[, mapping:=bsInferredMapping[regime]]
    setkey(bsInferredTreeHDDtNodes, endNode)

    bsInferredMeans <- try(PCMMean(
      bsInferredTreeHD, bsInferredModel, internal = TRUE), silent = TRUE)
    bsInferredVar <- try(PCMVar(
      bsInferredTreeHD, bsInferredModel, internal = TRUE)$Wii,
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

      attr(fit$inferredModel, "tree") <- NULL

      if(!keepBootstrapData) {
        attr(fit$inferredModel, "X") <- attr(fit$inferredModel, "SE") <- NULL
      }

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

