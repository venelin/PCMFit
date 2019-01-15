#' @export
GuessInitVecParams <- function(
  o, k = PCMNumTraits(o), R = PCMNumRegimes(o), n = 1L,
  argsPCMParamLowerLimit = NULL,
  argsPCMParamUpperLimit = NULL,
  X = NULL,
  tree = NULL,
  SE = NULL,
  tableAnc = NULL,
  returnWithinBoundsOnly = TRUE,
  ...) {
  UseMethod("GuessInitVecParams", o)
}

#' @export
GuessInitVecParams.PCM <- function(
  o, k = PCMNumTraits(o), R = PCMNumRegimes(o), n = 1L,
  argsPCMParamLowerLimit = NULL,
  argsPCMParamUpperLimit = NULL,
  X = NULL,
  tree = NULL,
  SE = NULL,
  tableAnc = NULL,
  returnWithinBoundsOnly = TRUE,
  ...) {

  res <- PCMParamRandomVecParams(
    o = o, k = k, R = R, n = n,
    argsPCMParamLowerLimit = argsPCMParamLowerLimit,
    argsPCMParamUpperLimit = argsPCMParamUpperLimit)

  if( !is.null(X) && "X0" %in% names(o) && !is.Fixed(o$X0) ) {
    vec.o <- PCMParamGetShortVector(o)
    vec.along.o <- seq_along(vec.o)

    o.along.o <- o

    PCMParamLoadOrStore(
      o.along.o, vec.along.o, offset = 0,
      k = PCMNumTraits(o), R = PCMNumRegimes(o), load = TRUE)

    if(is.null(tree)) {
      # no tree supplied, so use the population grand mean to set X0 in the
      # random parameter vectors
      grandMean <- rowMeans(X, na.rm = TRUE)
      grandMeanMask <- grandMean + 8108.20
      o.along.o$X0[] <- grandMeanMask
      v.along.o <- PCMParamGetShortVector(o.along.o)
      idxX0 <- match(grandMeanMask, v.along.o)
      idxX0 <- idxX0[!is.na(idxX0)]
      v.along.o[idxX0] <- v.along.o[idxX0] - 8108.20
      o.along.o$X0[] <- grandMean

      res[, idxX0] <- matrix(
        v.along.o[idxX0], nrow = n, ncol = length(idxX0), byrow = TRUE)
    } else {
      # a tree has been supplied, so calculate a weighted mean and assign it to
      # X0 in the random parameter vectors
      if(is.null(tableAnc)) {
        tableAnc <- PCMTreeTableAncestors(tree)
      }
      rootDists <- PCMTreeNodeTimes(tree)

      N <- PCMTreeNumTips(tree)

      daughtersRoot <- PCMTreeGetDaughters(tree, N + 1L)

      weightsDaughtersRoot <- 1 / rootDists[daughtersRoot]
      weightsDaughtersRoot <- weightsDaughtersRoot / sum(weightsDaughtersRoot)

      tipsFromDaughtersRoot <- lapply(daughtersRoot, function(j) {
        which(tableAnc[seq_len(PCMTreeNumTips(tree)), j] > 0)
      })

      regimesDaughtersRoot <- PCMTreeGetRegimesForNodes(tree, daughtersRoot)

      names(regimesDaughtersRoot) <- names(tipsFromDaughtersRoot) <-
        as.character(daughtersRoot)

      tipsFromDaughtersRootInSameRegime <-
        lapply(as.character(daughtersRoot), function(j) {
          intersect(tipsFromDaughtersRoot[[j]],
                    PCMTreeGetTipsInRegime(tree, regimesDaughtersRoot[[j]]))
        })

      names(tipsFromDaughtersRootInSameRegime) <- as.character(daughtersRoot)

      meansSubtreesDaughtersRoot <- matrix(
        sapply(daughtersRoot, function(j) {
          if(j <= N) {
            # j is a tip so use its value as is
            X[, j]
          } else {
            # j is an internal node. Calculate a weighted mean of the values
            # for the tips descending from j and having the same regime as j.
            # The weight of tip is is proportional to
            # 1/(rootDists[i] - rootDists[j]), that is, tips closer to j have
            # higher weight in calculating the mean.

            tips <- tipsFromDaughtersRootInSameRegime[[as.character(j)]]
            weightsTips <- 1 / (rootDists[tips] - rootDists[j])
            weightsTips <- weightsTips / sum(weightsTips)

            colSums(weightsTips * t(X[, tips, drop=FALSE]), na.rm = TRUE)
          }
        }),
        nrow = nrow(X),
        ncol = length(daughtersRoot))

      weightedMean <- colSums(
        weightsDaughtersRoot * t(meansSubtreesDaughtersRoot),
        na.rm = TRUE)

      weightedMeanMask <- weightedMean + 8108.20

      o.along.o$X0[] <- weightedMeanMask
      v.along.o <- PCMParamGetShortVector(o.along.o)
      idxX0 <- match(weightedMeanMask, v.along.o)
      idxX0 <- idxX0[!is.na(idxX0)]
      v.along.o[idxX0] <- v.along.o[idxX0] - 8108.20
      o.along.o$X0[] <- weightedMean

      res[, idxX0] <- matrix(
        v.along.o[idxX0], nrow = n, ncol = length(idxX0), byrow = TRUE)
    }

  }
  res
}


#' @export
GuessInitVecParams.OU <- function(
  o, k = PCMNumTraits(o), R = PCMNumRegimes(o), n = 1L,
  argsPCMParamLowerLimit = NULL,
  argsPCMParamUpperLimit = NULL,
  X = NULL,
  tree = NULL,
  SE = NULL,
  tableAnc = NULL,
  returnWithinBoundsOnly = TRUE,
  ...) {

  res <- NextMethod()

  listArgsRes <- c(as.list(environment()), list(...))
  do.call(GuessInitVecParamsOUInternal, listArgsRes)
}

#' @export
GuessInitVecParams.DOU <- function(
  o, k = PCMNumTraits(o), R = PCMNumRegimes(o), n = 1L,
  argsPCMParamLowerLimit = NULL,
  argsPCMParamUpperLimit = NULL,
  X = NULL,
  tree = NULL,
  SE = NULL,
  tableAnc = NULL,
  returnWithinBoundsOnly = TRUE,
  ...) {

  res <- NextMethod()

  listArgsRes <- c(as.list(environment()), list(...))
  do.call(GuessInitVecParamsOUInternal, listArgsRes)
}

GuessInitVecParamsOUInternal <- function(
  o, k = PCMNumTraits(o), R = PCMNumRegimes(o), n = 1L,
  argsPCMParamLowerLimit = NULL,
  argsPCMParamUpperLimit = NULL,
  X = NULL,
  tree = NULL,
  SE = NULL,
  tableAnc = NULL,
  varyTheta = FALSE,
  returnWithinBoundsOnly = TRUE,
  res,
  ...) {

  if(!is.null(X) && !is.null(tree)) {
    vec.o <- PCMParamGetShortVector(o)
    sd.along.o <- vec.along.o <- seq_along(vec.o)

    o.along.o <- o

    PCMParamLoadOrStore(
      o.along.o, vec.along.o, offset = 0,
      k = PCMNumTraits(o), R = PCMNumRegimes(o), load = TRUE)

    regimes <- as.character(PCMRegimes(o))
    rootDists <- PCMTreeNodeTimes(tree)
    idxsThetas <- c()

    for(r in seq_len(PCMNumRegimes(o))) {
      reg <- regimes[r]

      tipsInRegime <- PCMTreeGetTipsInRegime(tree, reg)
      weightsTipsInRegime <- rootDists[tipsInRegime]
      weightsTipsInRegime <- weightsTipsInRegime / sum(weightsTipsInRegime)

      meanTipsInRegime <- colSums(
        weightsTipsInRegime * t(X[, tipsInRegime, drop=FALSE]),
        na.rm = TRUE)
      varTipsInRegime <- colSums(
        weightsTipsInRegime * t(X[, tipsInRegime, drop=FALSE] - meanTipsInRegime)^2,
        na.rm = TRUE)

      maskNum <- 8108.20 + rnorm(1)

      meanTipsInRegimeMask <- meanTipsInRegime + maskNum

      o.along.o$Theta[, reg] <- meanTipsInRegimeMask
      v.along.o <- PCMParamGetShortVector(o.along.o)
      idxTheta <- match(meanTipsInRegimeMask, v.along.o)
      idxTheta <- idxTheta[!is.na(idxTheta)]
      vec.along.o[idxTheta] <- v.along.o[idxTheta] - maskNum


      o.along.o$Theta[, reg] <- sqrt(varTipsInRegime)
      s.along.o <- PCMParamGetShortVector(o.along.o)
      sd.along.o[idxTheta] <- s.along.o[idxTheta]

      idxsThetas <- c(idxsThetas, idxTheta)
    }

    if(varyTheta) {
      res[, idxsThetas] <- do.call(
        cbind,
        lapply(idxsThetas, function(i) {
          rnorm(n, vec.along.o[i], sd.along.o[i])
        }))
    } else {
      res[, idxsThetas] <- matrix(
        vec.along.o[idxsThetas],
        nrow = n,
        ncol = length(idxsThetas),
        byrow = TRUE)
    }

    if(returnWithinBoundsOnly) {
      # need to remove the parameters that go out of the lower-upper bound
      lowerModel <- do.call(PCMParamLowerLimit, c(list(o), argsPCMParamLowerLimit))
      lowerVecParams <- PCMParamGetShortVector(lowerModel)

      upperModel <- do.call(PCMParamUpperLimit, c(list(o), argsPCMParamUpperLimit))
      upperVecParams <- PCMParamGetShortVector(upperModel)

      withinBounds <- as.logical(
        sapply(seq_len(nrow(res)), function(i) {
          isTRUE(all(res[i,] >= lowerVecParams)) &&
            isTRUE(all(res[i,] <= upperVecParams))
        }))
        res <- res[withinBounds, , drop = FALSE]
    }
  }

  res
}

#' @export
GuessInitVecParams.MixedGaussian <- function(
  o, k = PCMNumTraits(o), R = PCMNumRegimes(o), n = 1L,
  argsPCMParamLowerLimit = NULL,
  argsPCMParamUpperLimit = NULL,
  X = NULL,
  tree = NULL,
  SE = NULL,
  tableAnc = NULL,
  returnWithinBoundsOnly = TRUE,
  ...) {

  res <- NextMethod()

  argsExtra <- list(...)
  if("varyTheta" %in% names(argsExtra)) {
    varyTheta <- argsExtra$varyTheta
  } else {
    varyTheta = FALSE
  }

  if(!is.null(X) && !is.null(tree)) {
    vec.o <- PCMParamGetShortVector(o)
    sd.along.o <- vec.along.o <- seq_along(vec.o)

    o.along.o <- o

    PCMParamLoadOrStore(
      o.along.o, vec.along.o, offset = 0,
      k = PCMNumTraits(o), R = PCMNumRegimes(o), load = TRUE)

    regimes <- as.character(PCMRegimes(o))
    rootDists <- PCMTreeNodeTimes(tree)
    idxsThetas <- c()

    for(r in seq_len(PCMNumRegimes(o))) {
      reg <- regimes[r]
      if(inherits(o[[reg]], "OU") || inherits(o[[reg]], "DOU")) {
        tipsInRegime <- PCMTreeGetTipsInRegime(tree, reg)

        weightsTipsInRegime <- rootDists[tipsInRegime]
        weightsTipsInRegime <- weightsTipsInRegime / sum(weightsTipsInRegime)

        meanTipsInRegime <- colSums(
          weightsTipsInRegime * t(X[, tipsInRegime, drop = FALSE]), na.rm = TRUE)
        varTipsInRegime <- colSums(
          weightsTipsInRegime * t(X[, tipsInRegime, drop = FALSE] - meanTipsInRegime)^2,
          na.rm = TRUE)

        maskNum <- 8108.20 + rnorm(1)

        meanTipsInRegimeMask <- meanTipsInRegime + maskNum

        o.along.o[[reg]]$Theta[] <- meanTipsInRegimeMask
        v.along.o <- PCMParamGetShortVector(o.along.o)
        idxTheta <- match(meanTipsInRegimeMask, v.along.o)
        idxTheta <- idxTheta[!is.na(idxTheta)]
        vec.along.o[idxTheta] <- v.along.o[idxTheta] - maskNum

        o.along.o[[reg]]$Theta[] <- sqrt(varTipsInRegime)
        s.along.o <- PCMParamGetShortVector(o.along.o)
        sd.along.o[idxTheta] <- s.along.o[idxTheta]

        idxsThetas <- c(idxsThetas, idxTheta)
      }
    }

    if(varyTheta) {
      res[, idxsThetas] <- do.call(
        cbind,
        lapply(idxsThetas, function(i) {
          rnorm(n, vec.along.o[i], sd.along.o[i])
        }))
    } else {
      res[, idxsThetas] <- matrix(
        vec.along.o[idxsThetas],
        nrow = n,
        ncol = length(idxsThetas),
        byrow = TRUE)
    }

    if(returnWithinBoundsOnly) {
      # need to remove the parameters that go out of the lower-upper bound
      lowerModel <- do.call(PCMParamLowerLimit, c(list(o), argsPCMParamLowerLimit))
      lowerVecParams <- PCMParamGetShortVector(lowerModel)

      upperModel <- do.call(PCMParamUpperLimit, c(list(o), argsPCMParamUpperLimit))
      upperVecParams <- PCMParamGetShortVector(upperModel)

      withinBounds <- sapply(1:nrow(res), function(i) {
        isTRUE(all(res[i,] >= lowerVecParams)) &&
          isTRUE(all(res[i,] <= upperVecParams))
      })
      res <- res[withinBounds, , drop = FALSE]
    }
  }
  res
}

