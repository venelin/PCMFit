#' Load a mixed-regime Gaussian model from a fit vector
#' @noRd
#' @keywords internal
PCMLoadMixedGaussianFromFitVector <- function(
  fitVector, modelTypes, k,
  remapModelTypeIndicesInFitVector = seq_along(modelTypes), ...) {
  # the last entries in fitVector are in the following order from left to right:
  # numNumericParams, logLik, df, nobs, score;
  # the first elements from 1 to numNumericParams are the actual numeric
  # parameters; the entries that follow are between the numeric parameters and
  # numNumericParams. These must be a pair number 2R, where R is the number of
  # regimes in the tree. The first R of these are the starting nodes for the R
  # regimes (starting from the root-node, which is always the starting node for
  # the first regime); The second R are the indices that map modelTypes to the
  # regimes on the tree.
  last <- length(fitVector)
  v_score <- fitVector[last]
  nobs <- fitVector[last - 1]
  df <- fitVector[last - 2]
  ll <- fitVector[last - 3]
  numNumericParams <- fitVector[last - 4]

  mappingModelsToRegimes <- matrix(
    as.integer(fitVector[(numNumericParams+1):(last-5)]), nrow = 2, byrow = TRUE)

  model <- MixedGaussian(
    k = k, modelTypes = modelTypes,
    mapping = remapModelTypeIndicesInFitVector[mappingModelsToRegimes[2, ]], ...)

  PCMParamLoadOrStore(model, fitVector, offset = 0, load = TRUE)

  attr(model, "startingNodesRegimes") <- mappingModelsToRegimes[1,]
  attr(model, "ll") <- ll
  attr(model, "df") <- df
  attr(model, "nobs") <- nobs
  attr(model, "score") <- v_score
  model
}

#' Match a model mapping vector against a vector of model types
#' @param modelMapping a character or integer vector
#' @param modelTypes a character vector with valid model-class-names.
#' @return a character vector with elements from modelTypes or stops with an error
#' @export
MatchModelMapping <- function(modelMapping, modelTypes) {
  if(is.character(modelMapping)) {
    m <- match(modelMapping, modelTypes)
    if( any(is.na(m)) ) {
      stop("MatchModelMapping:: some of the models in modelMapping could not be matched against model-types in tableFits (", toString(modelMapping[which(is.na(m))]), ")")
    }
  } else if( is.integer(modelMapping) ) {
    m <- match(modelMapping, seq_along(modelTypes))
    if( any(is.na(m)) ) {
      stop("MatchModelMapping:: some of the integer models in modelMapping could not be matched against model-indices in tableFits (", toString(modelMapping[which(is.na(m))]), ")")
    }
  } else {
    stop(
      "MatchModelMapping:: modelMapping should be character or integer (",
      toString(modelMapping), ")")
  }
  unname(modelTypes[m])
}


#' Compose a MixedGaussian model
#'
#' @keywords internal
#' @param tree,startingNodesRegimes,modelTypes,k,R,mapping,argsMixedGaussian,tableFits,modelTypesInTableFits,tableAncestors,verbose Internal use.
#'
#' @details This is an internal function that needs to be exported solely for
#' technical reasons. Not intended to be called by the end user.
#'
#' @return a PCM object.
#' @importFrom PCMBase PCMTreeExtractClade is.Fixed is.PCM PCMParamLoadOrStore PCMParamGetShortVector
#' @export
ComposeMixedGaussianFromFits <- function(
  tree, startingNodesRegimes, modelTypes, k, R, mapping,
  argsMixedGaussian,
  tableFits, modelTypesInTableFits,
  tableAncestors = NULL, verbose = FALSE) {

  if(verbose) {
    cat("Composing a MixedGaussian model:\n")
    cat("startingNodesRegimes=c(", toString(startingNodesRegimes), ")\n")
    cat("mapping=c(", toString(mapping), ")\n")
  }

  # tableFits can be supplied from a previous fit with different modelTypes which partially
  # overlap with modelTypes.
  # In this case, we need to remap the model indices in fitVector to the new modelTypes.
  remapModelTypeIndicesInFitVector <- match(modelTypesInTableFits, modelTypes)
  if(verbose) {
    cat("remapModelTypeIndicesInFitVector=c(",
        toString(remapModelTypeIndicesInFitVector), ");")
    cat("NAs indicate that some of modelTypesInTableFits are not present in modelTypes.\n")
  }

  # create a MixedGaussian model
  model <- do.call(
    MixedGaussian,
    c(list(k = k, modelTypes = modelTypes, mapping = mapping), argsMixedGaussian))

  # load listParInitOptim with the parameter vector from the ML fits to all clades
  # starting at startingNodesRegimes
  spec <- attr(model, "spec", exact = TRUE)
  subModelsLoaded <- rep(FALSE, length(mapping))

  for(r in seq_len(length(startingNodesRegimes))) {
    nr <- startingNodesRegimes[r]
    mr <- modelTypes[mapping[r]]
    tree_nr <- PCMTreeExtractClade(tree, nr, tableAncestors = tableAncestors)
    PCMTreeSetPartition(tree_nr)
    # model_mr <- do.call(
    #   MixedGaussian,
    #   c(list(k = k, modelTypes = modelTypes, mapping = mapping[r]), argsMixedGaussian))

    fit <- LookupFit(
      tree = tree_nr, modelTypes = modelTypes, modelMapping = mr,
      tableFits = tableFits)
    #fit <- LookupFit2(tree = tree_nr, model = model_nr, tableFits = tableFits)

    if(nrow(fit) == 0) {
      stop("ComposeMixedGaussianFromCladeFits:: no entry in tableFits for the given tree and modelMapping.")
    }

    fitVector <- fit$fitVector[[1]]

    if(verbose) {
      cat("Loading cladeFit for nr=", nr,"; mr=", mr, "\n")
    }

    fitModel <-
      do.call(PCMLoadMixedGaussianFromFitVector,
              c(list(fitVector = fitVector,
                     modelTypes = modelTypes,
                     k = k,
                     remapModelTypeIndicesInFitVector = remapModelTypeIndicesInFitVector),
                argsMixedGaussian))

    for(name in names(fitModel)) {
      if(verbose) {
        cat("Setting model or global parameter with name=", name, "\n")
      }

      if( name %in% names(model) &&
          is.Global(model[[name]]) &&
          is.Global(fitModel[[name]]) &&
          !is.Fixed(model[[name]]) ) {

        # set all global non-fixed parameters to the mean of their best fits
        # for the clades
        vecCurrent <- PCMParamGetShortVector(model[[name]], k = k, R = 1)
        vecFitModel <- PCMParamGetShortVector(fitModel[[name]], k = k, R = 1)

        if(! length(vecCurrent) == length(vecFitModel) ) {
          stop(
            paste0(
              "ComposeMixedGaussianFromCladeFits:: a global parameter ",
              name,
              " in a fit-model has a different short-vector length from the to-be fit model length(vecToFit)=",
              length(vecCurrent), ", length(vecFitModel)=", length(vecFitModel),
              "; class model to fit = ", mr,
              "; class fitted model = ", toString(class(fitModel[[name]]))))
        }
        vecCurrent <- vecCurrent + vecFitModel/R
        PCMParamLoadOrStore(model[[name]], vecCurrent, 0, k, 1, load = TRUE)

      } else if(is.PCM(spec[[name]]) && class(fitModel[[name]])[1] == mr) {
        if(!subModelsLoaded[r]) {
          model[[as.character(r)]] <- fitModel[[name]]
          subModelsLoaded[r] <- TRUE
        } else {
          stop(
            "ComposeMixedGaussianFromCladeFits:: submodel for regime ",
            r, " was already loaded from a best clade fit.")
        }
      } else {
        stop(
          "ComposeMixedGaussianFromCladeFits:: Found a member (",
          name, ") in fitModel starting from node (", nr, ") and with class '",
          class(fitModel[[name]]),
          "' which is neither a global parameter nor a model of the needed type (",
          mr, ").")
      }
    }
  }
  model
}

SaveCurrentResults <- function(listResults, filePrefix) {
  status <- try(
    save(listResults, file = paste0("Current_", filePrefix, ".RData")),
    silent = TRUE)
  if(inherits(status, 'try-error')) {
    warning(paste0("An error occurred while saving tableFits to file ",
                   paste0("Current_", filePrefix, ".RData"), " :", status))
  }
}

# combine fits from parallel tasks into a data.table and saves this data.table to
# an .RData file
CombineTaskResults <- function(..., envNCalls) {
  envNCalls$ncalls <- envNCalls$ncalls + 1
  data <- rbindlist(
    lapply(
      list(...),
      function(dt) {
        if(is.data.table(dt)) {
          dt
        } else {
          cat("ERROR in foreach task: no data.table returned: \n", toString(dt),
              "\n")
          NULL
        }
      }),
    use.names = TRUE)
  data
}

CleanTemporaryFitFiles <- function(filePrefix) {
  status <- try(
    {
      files <- list.files(".", pattern=paste0("^", filePrefix, ".*.RData"))
      if(length(files) > 0) {
        do.call(file.remove, as.list(files))
      }

    }, silent = TRUE)
  if(class(status) == 'try-error') {
    warning(paste0("An error occurred while cleaning up temporary data-files: ", status))
  }
}


