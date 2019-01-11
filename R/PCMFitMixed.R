#' Fit regime-assignments to (sub-)trees in a tree with different assigned model
#' types to each regime.
#'
#' @description This function performs multiple model fits of mixed regime models
#' (MixedGaussian) mapping different model-types (e.g. BM and OU) to different
#' regimes (colors) in a tree and testing different regime assignments to the
#' branches in the tree.
#' @importFrom foreach foreach when %do% %dopar% %:%
#' @importFrom data.table data.table rbindlist is.data.table setkey :=
#' @importFrom PCMBase PCMTreeSetLabels PCMTreeSetDefaultRegime PCMTreeEvalNestedEDxOnTree PCMTreeNumTips PCMTreeListCladePartitions PCMTreeListAllPartitions PCMTreeToString MixedGaussian PCMOptions PCMTreeTableAncestors PCMTreeSplitAtNode PCMTreeSetRegimes PCMGetVecParamsRegimesAndModels
#' @importFrom stats logLik coef AIC
#' @return an S3 object of class PCMFitModelMappings.
#'
#' @export
PCMFitMixed <- function(
  X, tree, modelTypes,
  SE = matrix(0.0, nrow(X), PCMTreeNumTips(tree)),

  generatePCMModelsFun = NULL,
  metaIFun = PCMInfo, positiveValueGuard = Inf,

  fitMappingsPrev = NULL,
  tableFitsPrev = fitMappingsPrev$tableFits,
  modelTypesInTableFitsPrev = NULL,

  listPartitions = NULL,
  minCladeSizes = 20L,

  maxCladePartitionLevel = 8L, maxNumNodesPerCladePartition = 1L,

  listAllowedModelTypesIndices = c("best-clade-2", "best-clade", "all"),

  scoreFun = AIC,

  argsMixedGaussian = NULL,

  argsConfigOptim1 = defaultArgsConfigOptim(numCallsOptim = 10),
  argsConfigOptim2 = defaultArgsConfigOptim(numCallsOptim = 4),

  listPCMOptions = PCMOptions(),

  skipFitWhenFoundInTableFits = TRUE,

  doParallel = FALSE,

  prefixFiles = "fits_",

  saveTempWorkerResults = TRUE,
  printFitVectorsToConsole = FALSE,
  setAttributes = TRUE,

  verbose = TRUE,
  debug = FALSE
) {

  if(is.list(listPartitions) || listPartitions == "all") {
    maxCladePartitionLevel = 1L
  }

  # Copy all arguments into a list
  # We establish arguments$<argument-name> as a convention for accessing the
  # original argument value.
  arguments <- as.list(environment())

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

    if(!is.list(arguments$listAllowedModelTypesIndices)) {
      # arguments$listAllowedModelTypesIndices is NULL or a character string
      listAllowedModelTypesIndices <-
        replicate(length(cladeRoots), seq_along(modelTypes), simplify = FALSE)
      names(listAllowedModelTypesIndices) <- as.character(cladeRoots)
    }

    if(verbose) {
      cat("Step 1: Performing fits on", length(cladeRoots), " clades; ",
          sum(sapply(listAllowedModelTypesIndices[as.character(cladeRoots)],
                     length)), " model mappings altogether...\n")
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
      listAllowedModelTypesIndices
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

    # update tableFits with the entries in fitsToClades which were not already
    # there.
    tableFits <- UpdateTableFits(tableFits, fitsToClades)

    SaveCurrentResults(list(tableFits = tableFits), filePrefix = prefixFiles)

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

  } else {
    stop("ERR:04131:PCMFit:PCMFitModelMappings.R:PCMFitModelMappings:: the tree is has fewer tips than than the min clade-size in minCladeSizes (",
         MIN_CLADE_SIZE,
         "). Try with smaller minCladeSizes.")
  }

  res <- list(
    options = PCMOptions(),
    tree = tree,
    X = X,
    SE = SE,
    hashCodeTree = resultStep2$hashCodeEntireTree,
    #tableFits = resultStep2$tableFits,
    tableFits = rbindlist(list(fitsToClades, resultStep2$fitsToTree)),
    queuePartitionRoots = resultStep2$queuePartitionRoots,
    mainLoopHistory = resultStep2$mainLoopHistory,
    arguments = arguments
  )
  class(res) <- "PCMFitModelMappings"
  res
}
