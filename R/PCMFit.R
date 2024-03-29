# Copyright 2018 Venelin Mitov
#
# This file is part of PCMFit
#
# PCMFit is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PCMFit is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PCMFit.  If not, see <http://www.gnu.org/licenses/>.

#' Fitting a PCM model
#' @inheritParams PCMBase::PCMLik
#' @param positiveValueGuard a real number (not necessarily positive) used during
#' the fit as a threshold for highly positive but likely incorrect log-likelihood
#' values. This argument is set to \code{Inf} by default and should be used only
#' as a last escape. The recommended way to of preventing such numerical errors
#' is by setting more stringent values for the runtime options
#' \code{PCMBase.Threshold.EV} and \code{PCMBase.Threshold.SV} (see
#' \code{\link{PCMOptions}}).
#' @param argsPCMParamLowerLimit,argsPCMParamUpperLimit named lists with
#' arguments passed to the functions \code{\link{PCMParamLowerLimit}} and
#' \code{\link{PCMParamUpperLimit}}, respectively. Default: \code{NULL}.
#' @param matParInit a matrix of any number of rows and p columns where, p is
#' the number of variable numerical parameters in the model
#' (equal to \code{PCMParamCount(model)}). Each row of this matrix specifies a
#' suggested starting location for the optim L-BFGS-B run. Default: NULL,
#' meaning that the initial parameters are to be chosen at random and/or using
#' calls to \code{\link{GuessInitVecParams}} function.
#' @param numRunifInitVecParams,numGuessInitVecParams integers specifying how
#' many parameter vectors should be drawn from a uniform distribution between
#' \code{PCMParamLowerLimit(model)} and \code{PCMParamUpplerLimit(model)}, and
#' how many parameter vectors should be generated by jittering the resulting
#' vector from a  call to \code{\link{GuessInitVecParams}}. Before starting the
#' optimization the model likelihood is evaluated at each of these vectors and
#' the top \code{numCallsOptim} vectors are chosen for starting locations for
#' optimization. The default settings are
#' \code{numRunifInitVecParams = if( is.null(matParInit) ) 1000L else 0L} and
#' \code{numGuessInitVecParams = if( is.null(matParInit) ) 100L else 0L}.
#' @param numCallsOptim integer specifying the maximum number of calls to
#' \code{\link{optim}} function. Default: 10. Note that this parameter would be
#' overwritten by a smaller \code{nrow(matParInit)} (if \code{matParInit} is
#' specified) or by a smaller number of generated initial parameter vectors that
#' satisfy the parameter limits (see also the arguments
#' \code{numRunifInitVecParams,numGuessInitVecParams} and
#' \code{argsPCMParamLowerLimit}, \code{argsPCMParamUpperLimit}).
#' @param control a list passed as control argument to \code{\link{optim}}.
#' Default: NULL.
#' @param doParallel logical indicating if optim calls should be executed in
#' parallel. Default: FALSE.
#' @param verbose logical indicating if information messages should be printed
#' to the console while running. Default: FALSE.
#' @return an object of class PCMFit
#' @importFrom PCMBase PCMCreateLikelihood PCMTree PCMInfo PCMParamCount PCMParamGetShortVector PCMParamLoadOrStore PCMParamLowerLimit PCMParamUpperLimit PCMParamRandomVecParams PCMOptions
#' @importFrom foreach foreach %do% %dopar%
#' @seealso \code{\link{PCMFitMixed}} \code{\link{PCMOptions}}
#'
#' @references
#' Mitov, V., Bartoszek, K., & Stadler, T. (2019). Automatic generation of
#'  evolutionary hypotheses using mixed Gaussian phylogenetic models.
#'  Proceedings of the National Academy of Sciences of the United States of
#'  America, 35, 201813823. http://doi.org/10.1073/pnas.1813823116
#'
#' Mitov, V., Bartoszek, K., Asimomitis, G., & Stadler, T. (2019). Fast
#'  likelihood calculation for multivariate Gaussian phylogenetic models with
#'  shifts. Theoretical Population Biology.
#'  http://doi.org/10.1016/j.tpb.2019.11.005

#' @export
PCMFit <- function(
  X, tree, model,
  SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  metaI = PCMInfo(X, tree, model, SE),
  positiveValueGuard = Inf,
  argsPCMParamLowerLimit = NULL,
  argsPCMParamUpperLimit = NULL,
  matParInit = NULL,
  numRunifInitVecParams = if( is.null(matParInit) ) 1000L else 0L,
  numGuessInitVecParams = if( is.null(matParInit) ) 100L else 0L,
  numCallsOptim = 10L,
  control = NULL,
  doParallel = FALSE,
  verbose = FALSE) {

  # Make sure tree is a PCMTree object (prevent potential errors in GuessInitVecParams)
  tree <- PCMTree(tree)

  lowerModel <- do.call(PCMParamLowerLimit, c(list(model), argsPCMParamLowerLimit))
  lowerVecParams <- PCMParamGetShortVector(lowerModel)

  upperModel <- do.call(PCMParamUpperLimit, c(list(model), argsPCMParamUpperLimit))
  upperVecParams <- PCMParamGetShortVector(upperModel)


  matParamsModel <- matrix(
    PCMParamGetShortVector(
      model, k = PCMNumTraits(model), R = PCMNumRegimes(model)), nrow = 1L)

  if(verbose) {
    cat("Generating", numRunifInitVecParams, "random init vectors...\n")
  }
  matParInitRunif <- PCMParamRandomVecParams(
    o = model,
    k = PCMNumTraits(model),
    R = PCMNumRegimes(model),
    n = numRunifInitVecParams,
    argsPCMParamLowerLimit = argsPCMParamLowerLimit,
    argsPCMParamUpperLimit = argsPCMParamUpperLimit
  )
  if(verbose) {
    cat("Generating", numGuessInitVecParams,
        "guessed init vectors with random jittering for ",
        numGuessInitVecParams - 1L, "of them...\n")
  }
  matParInitGuess <- GuessInitVecParams(
    o = model,
    k = PCMNumTraits(model),
    R = PCMNumRegimes(model),
    n = numGuessInitVecParams,
    argsPCMParamLowerLimit = argsPCMParamLowerLimit,
    argsPCMParamUpperLimit = argsPCMParamUpperLimit,
    X = X, tree = tree, SE = SE, varyParams = TRUE)

  matParInit <- rbind(
    matParInit,
    matParamsModel,
    matParInitRunif,
    matParInitGuess)

  if(verbose) {
    cat("Enforcing boundaries on init parameter vectors...\n")
  }
  EnforceBounds(matParInit, lowerVecParams, upperVecParams)

  if(nrow(matParInit) == 0L) {
    stop(
      paste(
        "In a call to PCMFit matParInit was empty. This can happen if you",
        "used only GuessInitVecParams to generate matParInit and all generated",
        "parameters were out of range. Consider using PCMParamRandomVecParams",
        "as well. "))
  }

  if(nrow(matParInit) > numCallsOptim) {
    if(verbose) {
      cat(
        "Evaluating likelihood at", nrow(matParInit), "parameter vectors...\n")
    }

    `%op%` <- if(isTRUE(doParallel) ||
                 (is.numeric(doParallel) && doParallel > 1)) `%dopar%` else `%do%`

    chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

    globalFuns <- unname(unlist(
      sapply(ls(.GlobalEnv), function(n) if(is.function(.GlobalEnv[[n]])) n else NULL)))

    # prevent 'no visible binding' warning
    is <- NULL

    valParInitOptim <- foreach(
      is = chunk(seq_len(nrow(matParInit)), 8L),
      .combine = c,
      .export = globalFuns) %op% {

        lik <- PCMCreateLikelihood(
          X = X, tree = tree, model = model, SE = SE,
          metaI = if(is.function(metaI)) {
            metaI(X = X, tree = tree, model = model, SE = SE)
          } else {
            metaI
          },
          positiveValueGuard = positiveValueGuard)

        sapply(is, function(i) lik(matParInit[i,]))
      }

    topVal <- order(valParInitOptim, decreasing = TRUE)[seq_len(numCallsOptim)]
    matParInit <- matParInit[topVal,, drop=FALSE]

    if(verbose) {
      cat(
        "Taking the top-", numCallsOptim,
        " parameter combinations sorted by decreasing log-likelihood value: ",
        toString(round(valParInitOptim[topVal], 2)), ".\n")
    }
  }

  res <- as.list(environment())
  # These objects tend to be very big. the lik function and the metaI object
  # can be recreated.
  res$lik <- res$chunk <- res$is <- res$`%op%` <- res$metaI <-
    res$globalFuns <- res$matParInit <-
    res$matParInitRunif <- res$matParInitGuess <-
    res$matParInitGuessVaryParams <- NULL

  res$PCMOptions <- PCMOptions()

  res <- c(
    res,
    RunOptim(
      X = X, tree = tree, model = model, SE = SE,
      metaI = metaI,
      positiveValueGuard = positiveValueGuard,
      parLower = lowerVecParams,
      parUpper = upperVecParams,
      matParInit = matParInit,
      control = control,
      doParallel = doParallel,
      verbose = verbose))

  if(!is.null(res$Optim)) {
    par <- res$Optim$par
    res$modelOptim <- model
    PCMParamLoadOrStore(res$modelOptim, vecParams = par, offset = 0, load = TRUE)
    res$logLikOptim <- res$Optim$value

  } else {
    res$modelOptim <- NULL
    res$logLikOptim <- as.double(PCMOptions()$PCMBase.Value.NA)
  }

  class(res) <- c("PCMFit", class(res))
  res
}

#' Retrieve the optimal PCM model from a PCMFit object
#' @param fit an object of S3 class \code{'PCMFit'}.
#' @seealso \code{\link{PCMFit}}.
#' @return a PCM model object.
#' @export
RetrieveBestModel <- function(fit) {
  if(is.PCMFit(fit)) {
    model <- fit$modelOptim
    attr(model, "tree") <- fit$tree
    attr(model, "X") <- fit$X
    attr(model, "SE") <- fit$SE
  } else {
    model <- NULL
  }
  model
}


#' Check if an object is a PCMFit.
#' @param object an R object.
#' @return logical.
#' @export
is.PCMFit <- function(object) {
  inherits(object, "PCMFit")
}

#' @importFrom PCMBase PCMOptions PCMTreeNumTips
#' @export
logLik.PCMFit <- function(object, ...) {
  if(!is.PCMFit(object)) {
    stop("object must inherit from class PCMFit.")
  }

  value <- as.double(PCMOptions()$PCMBase.Value.NA)

  if(!is.null(object$modelOptim)) {
    value <- object$logLikOptim
  }

  attr(value, "df") <- PCMParamCount(object$model,
                                     countRegimeChanges = TRUE,
                                     countModelTypes = TRUE)
  attr(value, "nobs") <- PCMTreeNumTips(object$tree)

  value
}

#' @importFrom PCMBase PCMGetVecParamsRegimesAndModels
#' @export
coef.PCMFit <- function(object, ...) {
  if(!is.PCMFit(object)) {
    stop("object must inherit from class PCMFit.")
  }

  model <- object$model

  if(!is.null(object$modelOptim)) {
    model <- object$modelOptim
  }

  par <- PCMGetVecParamsRegimesAndModels(model, object$tree, ...)
  c(par, numParam = PCMParamCount(model))
}

