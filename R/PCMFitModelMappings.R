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

  # copy all arguments into a list
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
    argumentsFitsToClades <- arguments[intersect(
      names(arguments), names(as.list(args(PCMFitModelMappingsToClades))))]

    argumentsFitsToClades$X <- X
    argumentsFitsToClades$tree <- tree
    argumentsFitsToClades$modelTypes <- modelTypes
    argumentsFitsToClades$SE <- SE
    argumentsFitsToClades$fitMappingsPrev <- NULL
    argumentsFitsToClades$tableFitsPrev <- tableFits
    argumentsFitsToClades$modelTypesInTableFitsPrev <- modelTypesInTableFitsPrev
    argumentsFitsToClades$minCladeSize <- MIN_CLADE_SIZE
    argumentsFitsToClades$argsConfigOptimAndMCMC <- argsConfigOptimAndMCMC1
    argumentsFitsToClades$preorderTree <- preorderTree
    argumentsFitsToClades$tableAncestors <- tableAncestors

    fitsToClades <- do.call(PCMFitModelMappingsToClades, argumentsFitsToClades)

    # update tableFits with the entries in fitsToClades which were not already there.
    tableFits <- UpdateTableFits(tableFits, fitsToClades)

    SaveCurrentResults(list(tableFits = tableFits), filePrefix = prefixFiles)

    # 2. Perform fits to clade-partitions with different model mappings
    # we need these variables throughout this step
    argumentsStep2 <- arguments[intersect(
      names(arguments), names(as.list(args(PCMFitRecursiveCladePartition))))]

    argumentsStep2$X <- X
    argumentsStep2$tree <- tree
    argumentsStep2$modelTypes <- modelTypes
    argumentsStep2$SE <- SE
    argumentsStep2$fitMappingsPrev <- NULL
    argumentsStep2$tableFitsPrev <- tableFits
    argumentsStep2$modelTypesInTableFitsPrev <- modelTypesInTableFitsPrev
    argumentsStep2$argsConfigOptimAndMCMC <- argsConfigOptimAndMCMC2
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
    tableFits = resultStep2$tableFits,
    queuePartitionRoots = resultStep2$queuePartitionRoots,
    mainLoopHistory = resultStep2$mainLoopHistory,
    arguments = arguments
  )
  class(res) <- "PCMFitModelMappings"
  res
}
