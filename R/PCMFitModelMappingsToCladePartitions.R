PCMFitModelMappingsToCladePartitions <- function(
  X, tree, modelTypes,
  SE = matrix(0.0, nrow(X), PCMTreeNumTips(tree)),

  listCladePartitions = NULL,
  listAllowedModelTypesIndices = NULL,

  generatePCMModelsFun = NULL,
  metaIFun = PCMInfo, positiveValueGuard = Inf,

  lik = NULL, prior = NULL, input.data = NULL, config = NULL,

  fitMappingsPrev = NULL,
  tableFitsPrev = fitMappingsPrev$tableFits,
  modelTypesInTableFitsPrev = NULL,

  skipFitWhenFoundInTableFits = TRUE,

  prefixFiles = "fits_",

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

  printFitVectorsToConsole = FALSE,

  doParallel = FALSE,

  verbose = TRUE,
  verbosePCMFit = FALSE,
  verboseComposeMixedGaussianFromFits = FALSE,
  verboseAdaptArgsConfigOptimAndMCMC = FALSE
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

  if(verbose) {
    cat("Fitting on ", length(listCladePartitions),
        " partitions of tree into clades using the model-fits on the clades as starting points...\n")
  }

  `%op%` <- if(doParallel) `%dopar%` else `%do%`

  envCombineTaskResults2 <- new.env()
  envCombineTaskResults2$ncalls <- 0

  fitsToTree <-
    foreach(
      cladePartition = listCladePartitions,
      iPartition = 1:length(listCladePartitions),
      .combine = function(...) rbindlist(list(...)),
      .multicombine = TRUE,
      .inorder = FALSE,
      .packages = (.packages()) ) %:%

    foreach(
      modelMapping = PCMIteratorMapping2(
        mapping = unlist(sapply(listAllowedModelTypesIndices[cladePartition], function(.) .[1])),
        modelTypes = seq_len(length(modelTypes)),
        allowedModelTypesIndices = lapply(listAllowedModelTypesIndices[cladePartition], sort)),

      .combine=function(...) {
        CombineTaskResults(
          ...,
          filePrefix = prefixFiles,
          envNCalls = envCombineTaskResults2)
      },
      .multicombine=TRUE,
      .inorder = FALSE,
      .errorhandling = "pass",
      .packages = (.packages()) ) %op% {
        try({
          # recreate tableAncestors under a different name, to avoid exporting a potentially big
          # tableAncestors
          if(!exists("tableAncestorsTree", .GlobalEnv)) {
            if(verbose) {
              cat("Creating tableAncestorsTree in .GlobalEnv during tree-fits\n")
            }
            assign("tableAncestorsTree", PCMTreeTableAncestors(tree, preorderTree), .GlobalEnv)
          }

          if(!exists("generatedPCMModels", .GlobalEnv) && !is.null(generatePCMModelsFun)) {
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

          PCMTreeSetRegimes(tree, nodes = as.integer(cladePartition))

          hashCodes <- HashCodes(
            tree = tree,
            modelTypes = modelTypes,
            startingNodesRegimesLabels = cladePartition,
            modelMapping = modelMapping)
          fit <- LookupFit(tableFits = tableFits, hashCodes = hashCodes)

          if(nrow(fit) == 1 && skipFitWhenFoundInTableFits) {
            dt.row <- fit
            if(verbose) {
              cat("Found a fit in tableFits on clade-partition: (",
                  toString(cladePartition),
                  "); mapping: (", toString(modelMapping), ")\n")
            }
            dt.row[, duplicated:=TRUE]
          } else {
            if(verbose) {
              cat("Performing ML fit on clade-partition: (",
                  toString(cladePartition),
                  "); mapping: (", toString(modelMapping), ")\n")
            }

            model <- ComposeMixedGaussianFromFits(
              tree = tree,
              startingNodesRegimes = as.integer(cladePartition),
              modelTypes = modelTypes,
              k = nrow(X),
              R = length(cladePartition),
              mapping = modelMapping,
              argsMixedGaussian = argsMixedGaussian,
              tableFits = tableFits,
              modelTypesInTableFits = modelTypesInTableFits,
              tableAncestors = tableAncestorsTree,
              verbose = verboseComposeMixedGaussianFromFits)

            fit <- PCMFit(
              X = X, tree = tree, model = model, SE = SE,
              metaI = metaIFun, positiveValueGuard = positiveValueGuard,
              lik = lik, prior = prior, input.data = input.data,
              config = config,
              argsPCMParamLowerLimit = argsPCMParamLowerLimit,
              argsPCMParamUpperLimit = argsPCMParamUpperLimit,
              argsPCMParamLoadOrStore = argsPCMParamLoadOrStore,
              argsConfigOptimAndMCMC = AdaptArgsConfigOptimAndMCMC(
                model,
                argsPCMParamLowerLimit = argsPCMParamLowerLimit,
                argsPCMParamUpperLimit = argsPCMParamUpperLimit,
                argsPCMParamLoadOrStore = argsPCMParamLoadOrStore,
                argsConfigOptimAndMCMC = argsConfigOptimAndMCMC,
                numJitterRootRegimeFit = numJitterRootRegimeFit,
                sdJitterRootRegimeFit = sdJitterRootRegimeFit,
                numJitterAllRegimeFits = numJitterAllRegimeFits,
                sdJitterAllRegimeFits = sdJitterAllRegimeFits,
                verbose = verboseAdaptArgsConfigOptimAndMCMC),
              verbose = verbosePCMFit)

            ll <- unname(logLik(fit))
            v_aic = unname(AIC(fit))
            vec <- c(
              coef(fit),
              logLik = ll,
              df = attr(ll, "df"),
              nobs = attr(ll, "nobs"),
              aic = v_aic)

            dt.row <- data.table(
              treeEDExpression = treeEDExpression,
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
              aic = v_aic,
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

          dt.row
        }, silent=FALSE)

      } # end of nested foreach body

  CleanTemporaryFitFiles(filePrefix = prefixFiles)
  fitsToTree
}
