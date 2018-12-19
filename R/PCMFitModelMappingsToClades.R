#' @importFrom foreach foreach when %do% %dopar% %:%
#' @importFrom data.table data.table rbindlist is.data.table setkey :=
#' @importFrom PCMBase PCMTreeSetLabels PCMTreeSetDefaultRegime PCMTreeEvalNestedEDxOnTree PCMTreeNumTips PCMTreeListCladePartitions PCMTreeToString MixedGaussian PCMOptions PCMTreeTableAncestors PCMTreeSplitAtNode PCMTreeSetRegimes PCMGetVecParamsRegimesAndModels
#' @importFrom stats logLik coef AIC
#' @export
PCMFitModelMappingsToClades <- function(
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

  minCladeSize = 25,

  listPCMOptions = PCMOptions(),

  argsMixedGaussian = NULL,
  argsPCMParamLowerLimit = NULL,
  argsPCMParamUpperLimit = NULL,
  argsPCMParamLoadOrStore = NULL,

  argsConfigOptimAndMCMC = NULL,

  preorderTree = NULL,
  tableAncestors = NULL,

  printFitVectorsToConsole = FALSE,

  doParallel = FALSE,

  verbose = TRUE,
  verbosePCMFit = FALSE
) {

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

  if(PCMTreeNumTips(tree) > minCladeSize) {

    if(is.null(preorderTree)) {
      preorderTree <- PCMTreePreorder(tree)
    }

    if(is.null(tableAncestors)) {
      tableAncestors <- PCMTreeTableAncestors(tree, preorder = preorderTree)
    }

    # 1. (fitsToClades) Perform a fit of each model-type to each clade
    # start from a list containing the trivial partition into one clade
    # equal to the whole tree
    cladeRoots <- c(PCMTreeNumTips(tree) + 1,
                    unlist(PCMTreeListCladePartitions(
                      tree = tree,
                      nNodes = 1,
                      minCladeSize = minCladeSize,
                      tableAncestors = tableAncestors)))
    if(verbose) {
      cat("Step 1: perform a fit on", length(cladeRoots), " clades x", length(modelTypes), " model-types...\n")
    }

    envCombineTaskResults1 <- new.env()
    envCombineTaskResults1$ncalls <- 0

    nodeLabelsTree <- PCMTreeGetLabels(tree)
    `%op%` <- if(doParallel) `%dopar%` else `%do%`

    fitsToClades <- foreach(
      clRoot = cladeRoots,
      clEDExpression =
        paste0("E(", treeEDExpression, ",", nodeLabelsTree[cladeRoots], ")"),
      .combine = function(...) rbindlist(list(...)),
      .multicombine = TRUE,
      .inorder = FALSE,
      .packages = (.packages()) ) %:%

      foreach(
        modelMapping = 1:length(modelTypes),
        .combine = function(...) {
          CombineTaskResults(
            ...,
            filePrefix = paste0(prefixFiles, "_clades_"),
            envNCalls = envCombineTaskResults1)
        },
        .multicombine = TRUE,
        .inorder = FALSE,
        .errorhandling = "pass",
        .packages = (.packages()) ) %op% {


          # recreate tableAncestors under a different name, to avoid
          # exporting a potentially big tableAncestors
          if(!exists("tableAncestorsTree", .GlobalEnv)) {
            if(verbose) {
              cat("Creating tableAncestorsTree in .GlobalEnv during clade-fits\n")
            }
            assign("tableAncestorsTree",
                   PCMTreeTableAncestors(tree, preorderTree), .GlobalEnv)
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

          # don't want the names in modelMapping
          modelMapping <- unname(modelMapping)

          options(listPCMOptions)

          k <- nrow(X)

          XSE <- rbind(X, SE)

          treeSplit <- PCMTreeSplitAtNode(
            tree, clRoot, tableAncestorsTree, XSE)
          clade <- treeSplit$clade
          Xclade <- treeSplit$Xclade[1:k, , drop = FALSE]
          SEclade <- treeSplit$Xclade[(k+1):(2*k), , drop = FALSE]

          PCMTreeSetDefaultRegime(clade, 1)

          hashCodes <- HashCodes(
            tree = clade,
            modelTypes = modelTypes,
            startingNodesRegimesLabels = PCMTreeGetLabels(clade)[
              PCMTreeNumTips(clade) + 1],
            modelMapping = modelMapping)

          fit <- LookupFit(tableFits = tableFits, hashCodes = hashCodes)

          if(nrow(fit) == 1 && skipFitWhenFoundInTableFits) {
            dt.row <- fit
            if(verbose) {
              cat("  Found a fit in tableFits for ",
                  modelTypes[modelMapping], " on clade starting at node ", clRoot, " ...\n")
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
            dt.row[, duplicated:=TRUE]

          } else {
            if(verbose) {
              cat("  Fitting ", modelTypes[modelMapping], " to clade starting at node ", clRoot, " ...\n")
            }

            # create an MixedGaussian model
            model <- do.call(MixedGaussian,
                             c(list(k = nrow(Xclade),
                                    modelTypes = modelTypes,
                                    mapping = modelMapping),
                               argsMixedGaussian))

            fit <- PCMFit(
              X = Xclade, tree = clade, model = model, SE = SEclade,
              metaI = metaIFun, positiveValueGuard = positiveValueGuard,

              lik = lik, prior = prior, input.data = input.data,
              config = config,
              argsPCMParamLowerLimit = argsPCMParamLowerLimit,
              argsPCMParamUpperLimit = argsPCMParamUpperLimit,
              argsPCMParamLoadOrStore = argsPCMParamLoadOrStore,
              argsConfigOptimAndMCMC = argsConfigOptimAndMCMC,
              verbose = verbosePCMFit)

            ll <- unname(logLik(fit))
            v_aic = unname(AIC(fit))
            vec <- c(coef(fit),
                     logLik = ll,
                     df = attr(ll, "df"),
                     nobs = attr(ll, "nobs"),
                     aic = v_aic)

            dt.row <- data.table(
              treeEDExpression = clEDExpression,
              hashCodeTree = hashCodes$hashCodeTree,
              hashCodeStartingNodesRegimesLabels =
                hashCodes$hashCodeStartingNodesRegimesLabels,

              hashCodeMapping = hashCodes$hashCodeMapping,
              # character names in global tree
              startingNodesRegimesLabels = list(
                PCMTreeGetLabels(clade)[PCMTreeNumTips(clade) + 1]),
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
                "\n",
                sep="")
          }
          dt.row
        }

    setkey(fitsToClades, hashCodeTree,hashCodeStartingNodesRegimesLabels,hashCodeMapping)

  }

  CleanTemporaryFitFiles(filePrefix = prefixFiles)

  fitsToClades
}
