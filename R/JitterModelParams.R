#' @importFrom PCMBase PCMParamGetShortVector PCMParamLoadOrStore
#' @importFrom stats rnorm
#' @export
JitterModelParams <- function(
  model,

  argsPCMParamLowerLimit = NULL,
  argsPCMParamUpperLimit = NULL,

  numJitterRootRegimeFit = 1000, sdJitterRootRegimeFit = 0.5,
  numJitterAllRegimeFits = 1000, sdJitterAllRegimeFits = 0.5,

  returnWithinBoundsOnly = TRUE,
  verbose = FALSE
) {

  matParams <-
    matrix(PCMParamGetShortVector(model), 1, PCMParamCount(model), byrow = TRUE)
  matParamsJitterRootCladeFit <- matParamsJitterAllCladeFits <- NULL

  # if there is more than one clade in the tree and numJitterRootRegimeFit > 0
  if( !is.null(model[["2"]]) && numJitterRootRegimeFit > 0 ) {
    vecParamIndex <- 1:ncol(matParams)
    modelIndexParams <- model
    PCMParamLoadOrStore(modelIndexParams, vecParamIndex, offset = 0, load = TRUE)
    vecParamIndexRootClade <- as.integer(PCMParamGetShortVector(modelIndexParams[["1"]]))
    matParamsJitterRootCladeFit <-
      matrix(
        matParams[1,],
        2 * numJitterRootRegimeFit,
        ncol(matParams),
        byrow=TRUE)
    for(j in vecParamIndexRootClade) {
      matParamsJitterRootCladeFit[, j] <-
        rnorm(
          2 * numJitterRootRegimeFit,
          mean = matParams[1, j],
          sd = sdJitterRootRegimeFit)
    }
  }

  if( numJitterAllRegimeFits > 0 ) {
    matParamsJitterAllCladeFits <-
      matrix(matParams[1, ],
             2 * numJitterAllRegimeFits,
             ncol(matParams),
             byrow=TRUE)
    for(j in 1:ncol(matParams)) {
      matParamsJitterAllCladeFits[, j] <- rnorm(2 * numJitterAllRegimeFits,
                                                mean = matParams[1, j],
                                                sd = sdJitterAllRegimeFits)
    }
  }

  if(!is.null(matParamsJitterRootCladeFit) || !is.null(matParamsJitterAllCladeFits)) {
    # need to remove the parameters that go out of the lower-upper bound
    lowerModel <- do.call(PCMParamLowerLimit, c(list(model), argsPCMParamLowerLimit))
    lowerVecParams <- PCMParamGetShortVector(lowerModel)

    upperModel <- do.call(PCMParamUpperLimit, c(list(model), argsPCMParamUpperLimit))
    upperVecParams <- PCMParamGetShortVector(upperModel)

    if( !is.null(matParamsJitterRootCladeFit) ) {
      if(returnWithinBoundsOnly) {
        withinBounds <- sapply(
          1:nrow(matParamsJitterRootCladeFit), function(i) {
            isTRUE(all(matParamsJitterRootCladeFit[i,] >= lowerVecParams)) &&
              isTRUE(all(matParamsJitterRootCladeFit[i,] <= upperVecParams))
          })
        matParamsJitterRootCladeFit <-
          matParamsJitterRootCladeFit[withinBounds, , drop = FALSE]
      }
      if(nrow(matParamsJitterRootCladeFit) > numJitterRootRegimeFit) {
        matParamsJitterRootCladeFit <-
          matParamsJitterRootCladeFit[1:numJitterRootRegimeFit, , drop = FALSE]
      }
      matParams <- rbind(matParams,
                                      matParamsJitterRootCladeFit)
    }

    if( !is.null(matParamsJitterAllCladeFits) ) {
      if(returnWithinBoundsOnly) {
        withinBounds <- sapply(
          1:nrow(matParamsJitterAllCladeFits), function(i) {
            isTRUE(all(matParamsJitterAllCladeFits[i,] >= lowerVecParams)) &&
              isTRUE(all(matParamsJitterAllCladeFits[i,] <= upperVecParams))
          })
        matParamsJitterAllCladeFits <-
          matParamsJitterAllCladeFits[withinBounds, , drop = FALSE]
      }
      if(nrow(matParamsJitterAllCladeFits) > numJitterAllRegimeFits) {
        matParamsJitterAllCladeFits <-
          matParamsJitterAllCladeFits[1:numJitterAllRegimeFits, , drop = FALSE]
      }
      matParams <- rbind(matParams,
                                      matParamsJitterAllCladeFits)
    }
  }

  matParams
}
