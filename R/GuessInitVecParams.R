#' A heuristic-based guess of the optimal parameters of a model
#'
#' @param o a PCM model object.
#' @param k integer denoting the number of traits in o.
#' Default: \code{PCMNumTraits(o)}.
#' @param n integer denoting the number of parameter vectors to generate.
#' @param argsPCMParamLowerLimit,argsPCMParamUpperLimit lists of parameters
#' passed to \code{PCMParamLowerLimit} and \code{PCMParamUpperLimit}
#' respectively. (see the argument \code{returnWithinBoundsOnly}).
#' @param tableAnc an ancestor table for tree. Default \code{NULL}.
#' @param varyParams logical indicating if fixed (FALSE) or varying (TRUE)
#' parameter values should be returned for each of the \code{n} vectors. The
#' default value is TRUE.
#' @param returnWithinBoundsOnly logical indication if the returned parameter
#' vectors should be within the lower and upper bound of the model.
#' @param res an n x p matrix where p is the number of parameters in o.
#' Internal use only.
#' @param ... additional arguments passed to implementations. Currently an
#' N x N matrix argument named treeVCVMat can be passed which is equal to the
#' phylogenetic covariance matrix for tree.
#' @return an n x p matrix where p is the number of parameters in o.
#'
#' @description This is an S3 generic function that returns a
#' \code{n}x\code{PCMParamCount(o)} matrix of parameter vectors.
#' Implementations try to deduce parameter values based on the passed tree
#' and data, that should be suitable for starting points for likelihood
#' optimization or MCMC inference.
#'
#' @examples
#' library(PCMBase)
#' library(PCMFit)
#' modelBM.ab <- modelBM.ab.Guess <- PCM(
#'   PCMDefaultModelTypes()["B"], k = 2, regimes = c("a", "b"))
#' modelBM.ab$X0[] <- c(5, 2)
#' # in regime 'a' the traits evolve according to two independent BM
#' # processes (starting from the global vecto X0).
#' modelBM.ab$Sigma_x[,, "a"] <- rbind(c(1.6, 0.01),
#'                                     c(0, 2.4))
#' # in regime 'b' there is a correlation between the traits
#' modelBM.ab$Sigma_x[,, "b"] <- rbind(c(2.4, 0.08),
#'                                     c(0.0, 0.8))
#' modelOU.ab.Guess <- PCM(
#'   PCMDefaultModelTypes()["F"], k = 2, regimes = c("a", "b"))
#'
#' modelMG.ab.Guess <- MixedGaussian(
#' k = 2, modelTypes = MGPMDefaultModelTypes(),
#' mapping = c(
#'   a = unname(MGPMDefaultModelTypes()["B"]),
#'   b = unname(MGPMDefaultModelTypes()["F"])))
#'
#'
#' param <- double(PCMParamCount(modelBM.ab))
#'
#' # load the current model parameters into param
#' PCMParamLoadOrStore(modelBM.ab, param, offset=0, load=FALSE)
#' print(param)
#'
#'
#' # make results reproducible
#' set.seed(2, kind = "Mersenne-Twister", normal.kind = "Inversion")
#'
#' # number of regimes
#' R <- 2
#'
#' # number of extant tips
#' N <- 100
#'
#' tree.a <- PCMTree(ape::rtree(n=N))
#' PCMTreeSetLabels(tree.a)
#' PCMTreeSetPartRegimes(
#'     tree.a, part.regime = c(`101` = "a"), setPartition = TRUE)
#'
#' lstDesc <- PCMTreeListDescendants(tree.a)
#' splitNode <- names(lstDesc)[which(sapply(lstDesc, length) > N/2 &
#'    sapply(lstDesc, length) < 2*N/3)][1]
#'
#' tree.ab <- PCMTreeInsertSingletons(
#'   tree.a, nodes = as.integer(splitNode),
#'   positions = PCMTreeGetBranchLength(tree.a, as.integer(splitNode))/2)
#' PCMTreeSetPartRegimes(
#'   tree.ab,
#'   part.regime = structure(
#'   c("a", "b"), names = as.character(c(N+1, splitNode))),
#'   setPartition = TRUE)
#'
#' traits <- PCMSim(tree.ab, modelBM.ab, modelBM.ab$X0)
#'
#' likFunBM <- PCMCreateLikelihood(traits, tree.ab, modelBM.ab.Guess)
#' likFunOU <- PCMCreateLikelihood(traits, tree.ab, modelOU.ab.Guess)
#' likFunMG <- PCMCreateLikelihood(traits, tree.ab, modelMG.ab.Guess)
#'
#' vecGuessBM <-
#'   GuessInitVecParams(modelBM.ab.Guess, X = traits, tree = tree.ab)
#' vecGuessOU <-
#'   GuessInitVecParams(modelOU.ab.Guess, X = traits, tree = tree.ab)
#' vecGuessMG <-
#'   GuessInitVecParams(modelMG.ab.Guess, X = traits, tree = tree.ab)
#'
#' # likelihood at true parameters used to generate the data
#' print(param)
#' likFunBM(param)
#'
#' # likelihood at guessed parameter values for the BM model
#' print(vecGuessBM)
#' likFunBM(vecGuessBM)
#'
#' # likelihood at guessed parameter values for the OU model
#' print(vecGuessOU)
#' likFunOU(vecGuessOU)
#'
#' # likelihood at guessed parameter values for the BM model
#' print(vecGuessMG)
#' likFunMG(vecGuessMG)
#'
#' # likelihood at completely random parameter values
#' vecRand <- PCMParamRandomVecParams(modelBM.ab.Guess, k = 2, R = 2)
#' likFunBM(vecRand)
#' @importFrom PCMBase PCMLik PCMNumRegimes PCMRegimes
#' @export
GuessInitVecParams <- function(
  o, regimes = as.character(PCMRegimes(o)), oParent = o, accessExpr = "",
  n = 1L,
  argsPCMParamLowerLimit = NULL,
  argsPCMParamUpperLimit = NULL,
  X = NULL,
  tree = NULL,
  SE = NULL,
  tableAnc = NULL,
  varyParams = FALSE,
  returnWithinBoundsOnly = TRUE,
  res = PCMParamRandomVecParams(
    o = oParent, k = PCMNumTraits(oParent), R = PCMNumRegimes(oParent), n = n,
    argsPCMParamLowerLimit = argsPCMParamLowerLimit,
    argsPCMParamUpperLimit = argsPCMParamUpperLimit),
  ...) {
  UseMethod("GuessInitVecParams", o)
}


#' @importFrom PCMBase PCMParamLocateInShortVector
#' @export
GuessInitVecParams.PCM <- function(
  o, regimes = as.character(PCMRegimes(o)), oParent = o, accessExpr = "",
  n = 1L,
  argsPCMParamLowerLimit = NULL,
  argsPCMParamUpperLimit = NULL,
  X = NULL,
  tree = NULL,
  SE = NULL,
  tableAnc = NULL,
  varyParams = FALSE,
  returnWithinBoundsOnly = TRUE,
  res = PCMParamRandomVecParams(
    o = oParent, k = PCMNumTraits(oParent), R = PCMNumRegimes(oParent), n = n,
    argsPCMParamLowerLimit = argsPCMParamLowerLimit,
    argsPCMParamUpperLimit = argsPCMParamUpperLimit),
  ...) {

  paramNames <- "X0"
  paramDefaultValues <- list()
  paramSDValues <- list()

  # Copy all arguments into a list
  # We establish arguments$<argument-name> as a convention for accessing the
  # original argument value.
  arguments <- c(as.list(environment()), list(...))

  res2 <- try(do.call(GuessParameters, arguments), silent = TRUE)

  if(!inherits(res2, "try-error")) {
    res <- res2
  }

  if(returnWithinBoundsOnly) {
    # need to remove the parameters that go out of the lower-upper bound
    lowerModel <- do.call(
      PCMParamLowerLimit, c(list(oParent), argsPCMParamLowerLimit))
    lowerVecParams <- PCMParamGetShortVector(lowerModel)

    upperModel <- do.call(
      PCMParamUpperLimit, c(list(oParent), argsPCMParamUpperLimit))
    upperVecParams <- PCMParamGetShortVector(upperModel)

    res <- EnforceBounds(res, lowerVecParams, upperVecParams)
  }

  res
}

#' @export
GuessInitVecParams.BM <- function(
  o, regimes = as.character(PCMRegimes(o)), oParent = o, accessExpr = "",
  n = 1L,
  argsPCMParamLowerLimit = NULL,
  argsPCMParamUpperLimit = NULL,
  X = NULL,
  tree = NULL,
  SE = NULL,
  tableAnc = NULL,
  varyParams = FALSE,
  returnWithinBoundsOnly = TRUE,
  res = PCMParamRandomVecParams(
    o = oParent, k = PCMNumTraits(oParent), R = PCMNumRegimes(oParent), n = n,
    argsPCMParamLowerLimit = argsPCMParamLowerLimit,
    argsPCMParamUpperLimit = argsPCMParamUpperLimit),
  ...) {

  paramNames <- c("X0", "Sigma_x")
  paramDefaultValues <- list()
  paramSDValues <- list()

  # Copy all arguments into a list
  # We establish arguments$<argument-name> as a convention for accessing the
  # original argument value.
  arguments <- c(as.list(environment()), list(...))

  res2 <- try(do.call(GuessParameters, arguments), silent = TRUE)

  if(!inherits(res2, "try-error")) {
    res <- res2
  }

  if(returnWithinBoundsOnly) {
    # need to remove the parameters that go out of the lower-upper bound
    lowerModel <- do.call(
      PCMParamLowerLimit, c(list(oParent), argsPCMParamLowerLimit))
    lowerVecParams <- PCMParamGetShortVector(lowerModel)

    upperModel <- do.call(
      PCMParamUpperLimit, c(list(oParent), argsPCMParamUpperLimit))
    upperVecParams <- PCMParamGetShortVector(upperModel)

    res <- EnforceBounds(res, lowerVecParams, upperVecParams)
  }

  res
}

#' @export
GuessInitVecParams.OU <- function(
  o, regimes = as.character(PCMRegimes(o)), oParent = o, accessExpr = "",
  n = 1L,
  argsPCMParamLowerLimit = NULL,
  argsPCMParamUpperLimit = NULL,
  X = NULL,
  tree = NULL,
  SE = NULL,
  tableAnc = NULL,
  varyParams = FALSE,
  returnWithinBoundsOnly = TRUE,
  res = PCMParamRandomVecParams(
    o = oParent, k = PCMNumTraits(oParent), R = PCMNumRegimes(oParent), n = n,
    argsPCMParamLowerLimit = argsPCMParamLowerLimit,
    argsPCMParamUpperLimit = argsPCMParamUpperLimit),
  ...) {

  paramNames <- c("X0", "H", "Theta", "Sigma_x", "Sigmae_x")
  paramDefaultValues <- list(H = 0.0)
  paramSDValues <- list(H = 0.02)

  # Copy all arguments into a list
  # We establish arguments$<argument-name> as a convention for accessing the
  # original argument value.
  arguments <- c(as.list(environment()), list(...))

  res2 <- try(do.call(GuessParameters, arguments), silent = TRUE)

  if(!inherits(res2, "try-error")) {
    res <- res2
  }

  if(returnWithinBoundsOnly) {
    # need to remove the parameters that go out of the lower-upper bound
    lowerModel <- do.call(
      PCMParamLowerLimit, c(list(oParent), argsPCMParamLowerLimit))
    lowerVecParams <- PCMParamGetShortVector(lowerModel)

    upperModel <- do.call(
      PCMParamUpperLimit, c(list(oParent), argsPCMParamUpperLimit))
    upperVecParams <- PCMParamGetShortVector(upperModel)

    res <- EnforceBounds(res, lowerVecParams, upperVecParams)
  }

  res
}

#' @export
GuessInitVecParams.DOU <- function(
  o, regimes = as.character(PCMRegimes(o)), oParent = o, accessExpr = "",
  n = 1L,
  argsPCMParamLowerLimit = NULL,
  argsPCMParamUpperLimit = NULL,
  X = NULL,
  tree = NULL,
  SE = NULL,
  tableAnc = NULL,
  varyParams = FALSE,
  returnWithinBoundsOnly = TRUE,
  res = PCMParamRandomVecParams(
    o = oParent, k = PCMNumTraits(oParent), R = PCMNumRegimes(oParent), n = n,
    argsPCMParamLowerLimit = argsPCMParamLowerLimit,
    argsPCMParamUpperLimit = argsPCMParamUpperLimit),
  ...) {

  paramNames <- c("X0", "H1", "H2", "Theta", "Sigmae_x", "Sigmae_x")
  paramDefaultValues <- list(H1 = 0.0, H2 = 0.0)
  paramSDValues <- list(H1 = 0.02, H2 = 0.02)

  # Copy all arguments into a list
  # We establish arguments$<argument-name> as a convention for accessing the
  # original argument value.
  arguments <- c(as.list(environment()), list(...))

  res2 <- try(do.call(GuessParameters, arguments), silent = TRUE)

  if(!inherits(res2, "try-error")) {
    res <- res2
  }

  if(returnWithinBoundsOnly) {
    # need to remove the parameters that go out of the lower-upper bound
    lowerModel <- do.call(
      PCMParamLowerLimit, c(list(oParent), argsPCMParamLowerLimit))
    lowerVecParams <- PCMParamGetShortVector(lowerModel)

    upperModel <- do.call(
      PCMParamUpperLimit, c(list(oParent), argsPCMParamUpperLimit))
    upperVecParams <- PCMParamGetShortVector(upperModel)

    res <- EnforceBounds(res, lowerVecParams, upperVecParams)
  }

  res
}


#' @importFrom PCMBase is.PCM is.Fixed
#' @export
GuessInitVecParams.MixedGaussian <- function(
  o, regimes = as.character(PCMRegimes(o)), oParent = o, accessExpr = "",
  n = 1L,
  argsPCMParamLowerLimit = NULL,
  argsPCMParamUpperLimit = NULL,
  X = NULL,
  tree = NULL,
  SE = NULL,
  tableAnc = NULL,
  varyParams = FALSE,
  returnWithinBoundsOnly = TRUE,
  res = PCMParamRandomVecParams(
    o = oParent, k = PCMNumTraits(oParent), R = PCMNumRegimes(oParent), n = n,
    argsPCMParamLowerLimit = argsPCMParamLowerLimit,
    argsPCMParamUpperLimit = argsPCMParamUpperLimit),
  ...) {

  # 1. first we guess all global parameters
  paramNames <- names(o)[!sapply(o, is.PCM) & !sapply(o, is.Fixed)]
  paramDefaultValues <- list(H = 0.0, H1 = 0.0, H2 = 0.0)
  paramSDValues <- list(H = 0.02, H1 = 0.02, H2 = 0.02)

  # Copy all arguments into a list
  # We establish arguments$<argument-name> as a convention for accessing the
  # original argument value.
  arguments <- c(as.list(environment()), list(...))

  res2 <- try(do.call(GuessParameters, arguments), silent = TRUE)

  if(!inherits(res2, "try-error")) {
    res <- res2
  }

  # 2. next, we call GuessInitVecParams for all submodels, which are not fixed
  for(reg in regimes) {
    if(!is.Fixed(o[[reg]])) {
      res2 <- try(
        GuessInitVecParams(
          o = o[[reg]], regimes = reg,
          oParent = o,
          accessExpr = paste0('[["', reg, '"]]'),
          n = n,
          argsPCMParamLowerLimit = argsPCMParamLowerLimit,
          argsPCMParamUpperLimit = argsPCMParamUpperLimit,
          X = X,
          tree = tree,
          SE = SE,
          tableAnc = tableAnc,
          varyParams = varyParams,
          returnWithinBoundsOnly = returnWithinBoundsOnly,
          res = res, ...),
        silent = TRUE)
      if(!inherits(res2, "try-error")) {
        res <- res2
      }
    }
  }

  if(returnWithinBoundsOnly) {
    # need to remove the parameters that go out of the lower-upper bound
    lowerModel <- do.call(
      PCMParamLowerLimit, c(list(oParent), argsPCMParamLowerLimit))
    lowerVecParams <- PCMParamGetShortVector(lowerModel)

    upperModel <- do.call(
      PCMParamUpperLimit, c(list(oParent), argsPCMParamUpperLimit))
    upperVecParams <- PCMParamGetShortVector(upperModel)

    res <- EnforceBounds(res, lowerVecParams, upperVecParams)
  }

  res
}

EnforceBounds <- function(vecs, lowerVecParams, upperVecParams) {
  hasNAs <- apply(vecs, 1, function(v) any(is.na(v)))
  vecs <- vecs[!hasNAs,, drop=FALSE]

  for(i in seq_len(nrow(vecs))) {
    tooSmall <- (vecs[i, ] < lowerVecParams)
    tooBig <- (vecs[i, ] > upperVecParams)
    vecs[i, tooSmall] <- lowerVecParams[tooSmall]
    vecs[i, tooBig] <- upperVecParams[tooBig]
  }
  vecs
}

#' @importFrom PCMBase PCMTreeGetDaughters PCMTreeGetTipsInPart
GuessX0 <- function(
  o, regimes,
  X = NULL,
  tree = NULL,
  SE = NULL,
  tableAnc = NULL,
  ...) {

  if( !is.null(X) && "X0" %in% names(o) && !is.Fixed(o$X0) ) {
    if(is.null(tree)) {
      # no tree supplied, so use the population grand mean to set X0 in the
      # random parameter vectors
      list(
        mean = apply(X, 1, mean, na.rm = TRUE),
        sd = apply(X, 1, sd, na.rm = TRUE))
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

      regimesDaughtersRoot <- PCMTreeGetPartsForNodes(tree, daughtersRoot)

      names(regimesDaughtersRoot) <- names(tipsFromDaughtersRoot) <-
        as.character(daughtersRoot)

      tipsFromDaughtersRootInSameRegime <-
        lapply(as.character(daughtersRoot), function(j) {
          intersect(tipsFromDaughtersRoot[[j]],
                    PCMTreeGetTipsInPart(tree, regimesDaughtersRoot[[j]]))
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

      list(
        mean = colSums(
          weightsDaughtersRoot * t(meansSubtreesDaughtersRoot),
          na.rm = TRUE),
        sd = apply(X, 1, sd, na.rm = TRUE))
    }
  } else {
    list(
      mean = o$X0[],
      sd = if(!is.null(o$X0)) {
        abs(o$X0[]) * 0.05
      } else {
        double(0L)
      })
  }
}

#' @importFrom PCMBase is.Global PCMTreeGetTipsInRegime PCM PCMDefaultModelTypes PCMVar PCMTreeGetTipsInRegime UpperChol PCMTreeVCV
GuessSigma_x <- function(
  o, regimes,
  X = NULL,
  tree = NULL,
  SE = NULL,
  tableAnc = NULL,
  ...) {

  Sigma_x <- SDSigma_x <- o$Sigma_x
  argsExtra <- list(...)

  if(!(is.null(Sigma_x) || is.Fixed(Sigma_x)) && !is.null(X) && !is.null(tree) ) {

    rootDists <- PCMTreeNodeTimes(tree, tipsOnly = TRUE)

    if(!is.null(argsExtra$treeVCVMat)) {
      treeVCVMat <- argsExtra$treeVCVMat
    } else {
      treeVCVMat <- PCMTreeVCV(tree)
    }

    CalculateRateMatrixBM <- function(X, C, rootDists) {
      N <- ncol(X)
      Cinv <- solve(C)
      weights <- rootDists / sum(rootDists)
      meanTipsInRegime <- colSums(weights * t(X), na.rm = TRUE)
      Y <- t(apply(t(X), 1, function(y) y - meanTipsInRegime))
      # a column matrix of 1s
      D <- matrix(1.0, nrow = N)
      # OLS solution
      BetaHat <- solve(t(D) %*% Cinv %*% D) %*% t(D) %*% Cinv %*% Y
      # rate matrix
      R <- t(Y - D%*%BetaHat) %*% Cinv %*% (Y - D%*%BetaHat) / N
      if(getOption("PCMBase.Transpose.Sigma_x", FALSE)) {
        chol(R)
      } else {
        UpperChol(R)
      }
    }

    SampleTips <- function(vecTips) {
      NMax <- getOption("PCMBase.MaxNForGuessSigma_x", 0.25)
      if(NMax <= 1) {
        # NMax is specified as a fraction of the number of tips
        sampSize <- NMax * length(vecTips)
        if(sampSize < 20 && length(vecTips) > 20) {
          # forcefully prevent using too few tips for the guess
          sampSize <- 20
        } else if(sampSize > 1000) {
          # forcefully prevent inverting to big matrices
          sampSize <- 1000
        }
        vecTips <- sample(vecTips, size = as.integer(sampSize))
      } else if(length(vecTips) > NMax) {
        vecTips <- sample(vecTips, size = NMax)
      }
    }

    if(is.Global(Sigma_x)) {
      # Sigma_x has a global scope for the entire tree
      tipsNA <- apply(X, 2, function(x) any(is.na(x)))
      idxTips <- seq_len(PCMTreeNumTips(tree))[!tipsNA]

      idxTips <- SampleTips(idxTips)

      Sigma_x[] <- CalculateRateMatrixBM(
        X[, idxTips, drop=FALSE],
        treeVCVMat[idxTips, idxTips, drop = FALSE],
        rootDists[idxTips])
      SDSigma_x[,] <- sqrt(abs(Sigma_x[,]) * 0.01)
    } else {
      # Sigma_x has a local scope for each regime in the model
      for(r in seq_len(PCMNumRegimes(o))) {
        reg <- regimes[r]

        tipsInRegime <- PCMTreeGetTipsInRegime(tree, reg)
        tipsNA <- apply(
          X[, tipsInRegime, drop=FALSE], 2, function(x) any(is.na(x)))
        tipsInRegime <- tipsInRegime[!tipsNA]

        tipsInRegime <- SampleTips(tipsInRegime)

        Sigma_x[,, r] <- CalculateRateMatrixBM(
          X[, tipsInRegime, drop=FALSE],
          treeVCVMat[tipsInRegime, tipsInRegime, drop = FALSE],
          rootDists[tipsInRegime])

        SDSigma_x[,, r] <- sqrt(abs(Sigma_x[,, r]) * 0.01)
      }
    }
  }

  list(mean = as.vector(Sigma_x), sd = as.vector(SDSigma_x))
}

GuessSigmae_x <- function(
  o, regimes,
  X = NULL,
  tree = NULL,
  SE = NULL,
  tableAnc = NULL,
  ...) {
  regimes <- as.character(PCMRegimes(o))

  Sigmae_x <- SDSigmae_x <- o$Sigmae_x

  if(!(is.null(Sigmae_x) || is.Fixed(Sigmae_x)) && !is.null(X)) {
    Sigmae_x[] <- 0.0
    SDSigmae_x[] <- 0.02

    if(is.Global(Sigmae_x)) {
      # Sigmae_x has a global scope for the entire tree
      for(j in seq_len(PCMNumTraits(o))) {
        Sigmae_x[j, j] <- sd(X[j,], na.rm = TRUE) * 0.025
        SDSigmae_x[j, j] <- sqrt(abs(Sigmae_x[j, j]) * 0.05)
      }
    } else {
      # Sigma_x has a local scope for each regime in the model
      for(r in seq_len(PCMNumRegimes(o))) {
        reg <- regimes[r]
        tipsInRegime <- PCMTreeGetTipsInRegime(tree, reg)

        for(j in seq_len(PCMNumTraits(o))) {
          Sigmae_x[j, j, r] <- sd(X[j, tipsInRegime], na.rm = TRUE) * 0.025
          SDSigmae_x[j, j, r] <- sqrt(abs(Sigmae_x[j, j, r]) * 0.05)
        }
      }
    }
  }

  list(mean = as.vector(Sigmae_x), sd = as.vector(SDSigmae_x))
}

#' @importFrom PCMBase PCMTreeNodeTimes
GuessTheta <- function(
  o, regimes,
  X = NULL,
  tree = NULL,
  SE = NULL,
  tableAnc = NULL,
  ...) {

  regimes <- as.character(PCMRegimes(o))

  Theta <- SDTheta <- o$Theta

  if(!(is.null(Theta) || is.Fixed(Theta)) && !is.null(tree) && !is.null(X)) {
    rootDists <- PCMTreeNodeTimes(tree, tipsOnly = TRUE)

    if(is.Global(Theta)) {
      weights <- rootDists / sum(rootDists)

      meanTips <- colSums(weights * t(X[, , drop=FALSE]), na.rm = TRUE)
      varTips <- colSums(weights * t(X[, , drop=FALSE] - meanTips)^2, na.rm = TRUE)

      Theta[] <- meanTips
      SDTheta[] <- sqrt(varTips)
    } else {
      # Sigma_x has a local scope for each regime in the model
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

        Theta[, r] <- meanTipsInRegime
        SDTheta[, r] <- sqrt(varTipsInRegime)
      }
    }
  }

  list(mean = as.vector(Theta), sd = as.vector(SDTheta))
}

GuessDefault <- function(
  name,
  value,
  sdValue,
  o, regimes,
  X = NULL,
  tree = NULL,
  SE = NULL,
  tableAnc = NULL,
  ...) {

  P <- SDP <- o[[name]]

  if(!(is.null(P) || is.Fixed(P)) && !is.null(X)) {
    P[] <- value
    SDP[] <- sdValue
  }

  list(mean = as.vector(P), sd = as.vector(SDP))
}

GuessParameters <- function(
  o, regimes, oParent = o, accessExpr = "",
  res,
  paramNames, paramDefaultValues = list(), paramSDValues = list(),
  n = 1L,
  argsPCMParamLowerLimit = NULL,
  argsPCMParamUpperLimit = NULL,
  X = NULL,
  tree = NULL,
  SE = NULL,
  tableAnc = NULL,
  varyParams = FALSE,
  returnWithinBoundsOnly = TRUE,
  ...) {

  for(name in paramNames) {
    if(name == "X0") {
      P <- GuessX0(o = o, regimes = regimes, X = X, tree = tree, SE = SE, tableAnc = tableAnc, ...)
    } else if(name == "Sigma_x") {
      P <- GuessSigma_x(o = o, regimes = regimes, X = X, tree = tree, SE = SE, tableAnc = tableAnc, ...)
    } else if(name == "Sigmae_x") {
      P <- GuessSigmae_x(o = o, regimes = regimes, X = X, tree = tree, SE = SE, tableAnc = tableAnc, ...)
    } else if(name == "Theta") {
      P <- GuessTheta(o = o, regimes = regimes, X = X, tree = tree, SE = SE, tableAnc = tableAnc, ...)
    } else {
      P <- GuessDefault(
        name, paramDefaultValues[[name]], paramSDValues[[name]],
        o = o, regimes = regimes, X = X, tree = tree, SE = SE, tableAnc = tableAnc, ...)
    }

    if(length(P$mean) > 0L) {
      i <- PCMParamLocateInShortVector(oParent, paste0(accessExpr, "$", name))
      iP <- which(!is.na(i))
      res[, iP] <- matrix(
        P$mean[i[iP]], nrow = n, ncol = length(iP), byrow = TRUE)
      if(varyParams && n > 1L) {
        res[-1L, iP] <- do.call(
          cbind, lapply(iP, function(j) rnorm(n-1L, P$mean[i[j]], P$sd[i[j]]) ) )
      }
    }
  }

  res
}
