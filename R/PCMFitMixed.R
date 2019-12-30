#' Search for an optimal mixed Gaussian phylogenetic model, given a tree, trait
#' measurements at its tips and a score function.
#'
#' @description A mixed Gaussian phylogenetic model (MGPM) represents a
#' Gaussian phylogenetic model with shifts in the underlying parameters and,
#' optionally, type of Gaussian stochastic process (e.g. shifts from a BM to an
#' OU model of evolution). The function \code{PCMFitMixed} implements a
#' recursive clade partition (RCP) search for an approximate information score
#' optimization. A deteailed description of the algorithm is provided in
#' Appendix A in Mitov et al. 2019a. For this documentation, it is important to
#' note that the algorithm proceeds in three steps as follows:
#' \describe{
#' \item{clade-fits}{Model type fits to all clades in the tree bigger than a specified
#' (see argument \code{minCladeSize}).}
#' \item{RCP}{Recursive clade partition search.}
#' \item{Round-robin}{Round-robin search for an optimal model type mapping in
#' the best partitions found during the RCP-step. This step is optional.}
#' }
#' @inheritParams PCMFit
#' @param doParallel logical indicating if the recursive clade partition search
#' should be executed in parallel. Default: \code{FALSE}.
#' @param verbose logical indicating if information messages should be printed
#' to the console while running. Default: \code{TRUE}.
#' @param modelTypes a named character string vector. Each such
#' string denotes a MGPM model class. The default setting is
#' \code{MGPMDefaultModelTypes()}, which corresponds to the model types
#' \code{A,B,C,D,E,F} as defined in Mitov et al. 2019a:
#' \describe{
#' \item{A. }{BM (H = 0, diagonal Σ): BM, uncorrelated traits.}
#' \item{B. }{BM (H = 0, symmetric Σ): BM, correlated traits.}
#' \item{C. }{OU (diagonal H, diagonal Σ): OU, uncorrelated traits.}
#' \item{D. }{OU (diagonal H, symmetric Σ): OU, correlated traits, but simple
#' (diagonal) selection strength matrix.}
#' \item{E. }{OU (symmetric H, symmetric Σ): An OU with nondiagonal symmetric H
#' and nondiagonal symmetric Σ.}
#' \item{F. }{OU (asymmetric H, symmetric Σ): An OU with nondiagonal asymmetric
#' H and nondiagonal symmetric Σ.}
#' }
#' See also the argument \code{subModels} and
#' \code{\link{MGPMDefaultModelTypes}}.
#' @param subModels a named character string vector with names and elements
#' among the names of \code{modelTypes}. This argument specifies a nesting
#' relationship between the model types: a model A is a sub-model of a
#' model B if A's parameters are a subset of B's parameters. During the RCP
#' search, the parameters of a sub-model fit to a clade in the tree are used as
#' a hint for the parameters of its encompassing models. The default setting is
#' \code{c(B = 'A', C = 'A', D = 'B', E = 'D', F = 'E')}. Note that this
#' setting is valid only for the default setting of the argument
#' \code{modelTypes}. Setting this argument to \code{NULL} will turn off the
#' use of nested models as hints.
#' @param argsMixedGaussian a list of arguments passed to
#' \code{\link{MixedGaussian}}. Default:
#' \code{Args_MixedGaussian_MGPMDefaultModelTypes()}.
#' @param generatePCMModelsFun A function generating PCM model types.
#' Specifying this function is only needed if the inference is to be done using #' custom candidate model types. Default setting:
#' \code{\link{PCMGenerateModelTypes}}. Note that this function is called in
#' each worker process and populates the global environment with
#' \code{\link{PCMParentClasses}} and \code{\link{PCMSpecify}} S3 methods. See
#' also \code{\link{PCMGenerateModelTypes}}.
#' @param metaIFun a metaI function needed to generate meta- and cache- objects
#' for the model objects generated during the RCP search. By default this is
#' set to \code{\link{PCMInfo}} from the package PCMBase. For faster execution,
#' it is recommended to use the function \code{PCMInfoCpp} from the package
#' PCMBaseCpp.
#' @param scoreFun a information score function such as AIC or BIC.
#' Default: \code{\link{AIC}}.
#' @param fitMappingsPrev an object of S3 class PCMFitModelMappings, returned
#' during a previous call to \code{PCMFitMixed}. Default: \code{NULL}. This can
#' be used to enrich a previous MGPM inference with new model types. Currently,
#' this feature is experimental. See also arguments \code{tableFitsPrev} and
#' \code{modelTypesInTableFitsPrev}.
#' @param tableFitsPrev a \code{\link{data.table}} object containing stored
#' MGPM fits from a previous (possibly unfinished) call to \code{PCMFitMixed}.
#' During a run \code{PCMFitMixed} stores all visited candidate MGPM models in
#' a file named 'Current_X.RData' where X is the value of the argument
#' \code{prefixFiles}. This file contains a list object of class
#' 'PCMFitModelMappings' named \code{'listResults'}. One of the elements in
#' \code{'listResults'} is a data.table \code{'tableFits'} that can be passed
#' as a value of this argument. Default: \code{fitMappingsPrev$tableFits}.
#' If a candidate model is found in tableFitsPrev, the previously inferred
#' parameters will be reused unless the argument
#' \code{skipFitWhenFoundInTableFits} is set to \code{FALSE}. See
#' also the argument \code{skipFitWhenFoundInTableFits}.
#' @param skipFitWhenFoundInTableFits logical. Default: \code{TRUE}.
#' @param modelTypesInTableFitsPrev model types used for the fits in
#' \code{'tableFitsPrev'}. Default: \code{NULL}, which means that the same
#' model types are used as specified in the argument \code{modelTypes}.
#' @param prefixFiles a character string used to name a 'Current_X.RData' file
#' stored during the inference, where X is to be replaced by the value of
#' \code{prefixFiles}. Note that the name
#' 'prefixFiles' can be misleading, because the actual file-name does not start
#' with it - this is due for compliance with existing code. Default:
#' \code{"fits_"}.
#' @param listPartitions,minCladeSizes,skipNodes arguments controlling the
#' search for an optimal shift point configuration. If not \code{NULL},
#' \code{listPartitions} should be specified as a list of integer vectors with
#' each such vector containing node indices in \code{tree} corresponding to a
#' shift-point configuration (partition). If \code{listPartitions} is
#' \code{NULL} (default) then the shift-point configuration are generated
#' according to the RCP search algorithm and the arguments \code{minCladeSizes}
#' and \code{skipNodes}. The argument \code{minCladeSizes} is an integer or a
#' vector of integers - if a single integer it specifies the the minimum number
#' of tips allowed in a part at al recursion levels of the RCP algorithm; if a
#' vector, then each element specifies the minimal part size at each level of
#' the algorithm (see also argument \code{maxCladePartitionLevel}). By default
#' this is set to \code{20}. The argument \code{skipNodes} is used in
#' combination with \code{minCladeSizes}. This is an integer or character
#' vector indicating the ids or labels of nodes that should not be used as
#' partition nodes. By default, this is an empty character vector. Default
#' setting for the three arguments:
#' \code{
#'   listPartitions = NULL,
#'   minCladeSizes = 20L,
#'   skipNodes = character()}. See also \code{maxCladePartitionLevel},
#'  \code{maxNumNodesPerCladePartition} and
#'  \code{listAllowedModelTypesIndices}.
#' @param maxCladePartitionLevel An integer limiting the number of recursive
#' levels of the RCP search algorithm. Default setting:
#' \code{if(is.null(listPartitions)) 8L else 1L}.
#' @param maxNumNodesPerCladePartition An integer controlling the number of
#' partition nodes in subtree-partitions during the RCP search. By default,
#' this is set to \code{1}, meaning that only single shifts within a parent part
#' will be searched. Setting this argument to a bigger value will make the
#' RCP search algorithm 'less greedy' but will cause a significant slow-down
#' in the search - many more partitions will be visited during the search.
#' @param listAllowedModelTypesIndices
#' Default setting: \code{c("best-clade-2", "best-clade", "all")}.
#' @param argsConfigOptim1,argsConfigOptim2,argsConfigOptim3 Arguments
#' controlling the calls to \code{\link{PCMFit}} during the three steps of the
#' RCP algorith.
#  Default setting:
#' \code{
#' argsConfigOptim1 = \link{DefaultArgsConfigOptim}(numCallsOptim = 10),
#' argsConfigOptim2 = \link{DefaultArgsConfigOptim}(numCallsOptim = 4),
#' argsConfigOptim3 = \link{DefaultArgsConfigOptim}(numCallsOptim = 10)
#' }.
#' @param maxNumRoundRobins number of round-robin iterations. By default this is
#' set ot 0, i.e. no round-robin step is performed.
#' @param maxNumPartitionsInRoundRobins maximum number (2 by default) of top
#' partitions tobe included in the round-robin step. This argument only has an
#' impact if \code{maxNumRoundRobins > 0}.
#' @param listPCMOptions a list of PCM runtime options to be specified prior
#' to running the fit. By default, this is set to \code{\link{PCMOptions}()}.
#' @param saveTempWorkerResults,printFitVectorsToConsole,debug
#' Debugging options controlling log-messages and storing of temporary results.
#' These are mostly for internal use. Default setting:
#' \code{saveTempWorkerResults, = TRUE, printFitVectorsToConsole = FALSE,
#' debug = FALSE}
#'
#' @importFrom foreach foreach when %do% %dopar% %:%
#' @importFrom data.table data.table rbindlist is.data.table setkey :=
#' @importFrom PCMBase PCMTree PCMTreeSetLabels PCMTreeSetPartition
#' PCMTreeEvalNestedEDxOnTree PCMTreeNumTips PCMTreeListCladePartitions
#' PCMTreeListAllPartitions PCMTreeToString MixedGaussian PCMOptions
#' PCMTreeTableAncestors PCMTreeSplitAtNode PCMGetVecParamsRegimesAndModels
#' MGPMDefaultModelTypes PCMGenerateModelTypes is.Transformable PCMTreeVCV
#' Args_MixedGaussian_MGPMDefaultModelTypes PCMTreeMatchLabels PCMTreePreorder
#' @importFrom stats logLik coef AIC
#'
#' @return an object of S3 class PCMFitModelMappings.
#' @seealso \code{\link{PCMFitMixed}} \code{\link{PCMOptions}}
#' @references
#' [Mitov et al. 2019a] Mitov, V., Bartoszek, K., & Stadler, T. (2019).
#' Automatic generation of evolutionary hypotheses using mixed Gaussian
#' phylogenetic models. Proceedings of the National Academy of Sciences of the
#' United States of America, 35, 201813823.
#' http://doi.org/10.1073/pnas.1813823116
#'
#' [Mitov et al. 2019b] Mitov, V., Bartoszek, K., Asimomitis, G., & Stadler, T. (2019). Fast
#'  likelihood calculation for multivariate Gaussian phylogenetic models with
#'  shifts. Theoretical Population Biology.
#'  http://doi.org/10.1016/j.tpb.2019.11.005
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
  skipNodes = character(),

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

  # Prevent no visible binding warnings:

  if( !is.null(listPartitions) ) {
    maxCladePartitionLevel = 1L
  }

  # Copy all arguments into a list
  # We establish arguments$<argument-name> as a convention for accessing the
  # original argument value.
  arguments <- as.list(environment())

  optionsBeforeCall <- options()[names(listPCMOptions)]

  do.call(options, listPCMOptions)


  tree <- PCMTree(tree)
  # this will set the nodelabels to the character representation of N+1, N+2, ..., M
  # the original node-labels are kept in arguments$tree.
  PCMTreeSetLabels(tree)
  PCMTreeSetPartition(tree)

  if(is.character(skipNodes) & length(skipNodes) > 0) {
    # get the converted node-labels to be skipped
    skipNodes <-
      PCMTreeGetLabels(tree)[PCMTreeMatchLabels(arguments$tree, skipNodes)]
  }

  colnames(X) <- dimnames(SE)[[length(dim(SE))]] <- as.character(seq_len(PCMTreeNumTips(tree)))

  tableFits <- InitTableFits(modelTypes,
                             fitMappingsPrev,
                             tableFitsPrev,
                             modelTypesInTableFitsPrev,
                             verbose = verbose)

  modelTypesInTableFits <- attr(tableFits, "modelTypes")

  arguments$tableFitsPrev <- arguments$fitMappingsPrev <- paste0(
    "To save space and avoid redundancy, tableFitsPrev and fitMappingsPrev\n",
    "are not stored in arguments. Most fits in tableFitsPrev are found\n",
    "in the member tableFits in the PCMFitModelMappings object. However some\n",
    "of these fits might have been replaced by equivalent model fits with a\n",
    "higher score encountered during this call to PCMFitMixed search.")

  if(length(minCladeSizes) < maxCladePartitionLevel) {
    minCladeSizes <- c(
      rep(as.integer(NA), maxCladePartitionLevel - length(minCladeSizes)),
      minCladeSizes)
  }
  MIN_CLADE_SIZE <- min(minCladeSizes, na.rm = TRUE)

  if(PCMTreeNumTips(tree) > MIN_CLADE_SIZE) {

    preorderTree <- PCMTreePreorder(tree)
    tableAncestors <- PCMTreeTableAncestors(tree, preorder = preorderTree)
    treeVCVMat <- PCMTreeVCV(tree)

    # 1. (fitsToClades) Perform a fit of each model-type to each clade
    if(is.null(arguments$listPartitions) || arguments$listPartitions == "all") {
      cladeRoots <- c(PCMTreeNumTips(tree) + 1,
                      unlist(PCMTreeListCladePartitions(
                        tree = tree,
                        nNodes = 1,
                        minCladeSize = MIN_CLADE_SIZE,
                        skipNodes = skipNodes,
                        tableAncestors = tableAncestors)))
    } else {
      cladeRoots = setdiff(
        unique(c(PCMTreeNumTips(tree) + 1,
                 unlist(arguments$listPartitions))),
        as.integer(skipNodes))
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
    argumentsFitsToClades$treeVCVMat <- treeVCVMat
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
    while(checkForBetterSubmodels &&
          checkForBetterSubmodelsIteration <= length(subModels)) {

      checkForBetterSubmodelsIteration <- checkForBetterSubmodelsIteration + 1L
      if(verbose) {
        cat(
          "Step 1.2, Iteration ", checkForBetterSubmodelsIteration,
          "(", Sys.time() ,"):",
          "Learning from sub-models, where the found max log-likelihood of a",
          "super-model was lower than the one of its sub-model...\n")
      }

      betterSubmodelFits <- UpdateCladeFitsUsingSubModels(
        cladeFits = fitsToClades,
        modelTypes = modelTypes,
        subModels = subModels,
        argsMixedGaussian = argsMixedGaussian,
        metaIFun = metaIFun,
        scoreFun = scoreFun,
        X = X, tree = tree, SE = SE,
        verbose = verbose)

      if(nrow(betterSubmodelFits$cladeFitsNew) > 0L) {

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
    argumentsStep2$treeVCVMat <- treeVCVMat
    argumentsStep2$SE <- SE
    argumentsStep2$skipNodes <- skipNodes
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

    # prevent 'no visible binding notes'
    score <- hashCodeStartingNodesRegimesLabels <-
      startingNodesRegimesLabels <- mapping <- df <- oldScore <- NULL

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

    resFitMappings <- list(
      arguments = arguments,
      options = listPCMOptions,
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

    SaveCurrentResults(resFitMappings, filePrefix = prefixFiles)

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

            resFitMappings <- list(
              arguments = arguments,
              options = listPCMOptions,
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
            SaveCurrentResults(resFitMappings, filePrefix = prefixFiles)
          }

          if(verbose) {
            cat(
              "> Step 3, iteration ", iRR, " (", Sys.time(), ")",
              "round robin loop for node position", pos,
              ", top candidate partitions/mappings after iterating over model types for this position:\n")
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
    stop("PCMFitMixed: the tree has fewer tips than the min clade-size in minCladeSizes (",
         MIN_CLADE_SIZE,
         "). Try with smaller minCladeSizes.")
  }

  resFitMappings <- list(
    arguments = arguments,
    options = listPCMOptions,
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

  resetOptions <- try(do.call(options, optionsBeforeCall), silent = TRUE)
  if(inherits(resetOptions, "try-error")) {
    warning(paste0(
      "PCMFitMixed:: Could not reset runtime options to the original values.",
      resetOptions))
  }

  resFitMappings
}
