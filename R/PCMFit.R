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


#' Fitting a PCM model to a given tree and data
#' @inheritParams PCMBase::PCMLik
#' @return an object of class PCMFit
#' @export
PCMFit <- function(X, tree, model,
                   SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
                   ...) {
  UseMethod("PCMFit", model)
}


#' @importFrom OptimMCMC runOptimAndMCMC configOptimAndMCMC
#' @importFrom PCMBase PCMCreateLikelihood PCMInfo PCMParamCount PCMParamGetShortVector PCMParamLoadOrStore PCMParamLowerLimit PCMParamUpperLimit PCMOptions
#' @export
PCMFit.PCM <- function(
  X, tree, model,
  SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  metaI = PCMInfo(X, tree, model, SE), positiveValueGuard = Inf,
  lik = NULL, prior = NULL, input.data = NULL, config = NULL,
  argsPCMParamLowerLimit = NULL,
  argsPCMParamUpperLimit = NULL,
  argsPCMParamLoadOrStore = NULL,
  argsConfigOptimAndMCMC = NULL,
  verbose = FALSE, ...) {

  if(is.null(lik)) {
    lik <- PCMCreateLikelihood(X, tree, model, SE, metaI, positiveValueGuard)
  }
  if(is.null(config)) {
    lowerModel <- do.call(PCMParamLowerLimit, c(list(model), argsPCMParamLowerLimit))
    lowerVecParams <- PCMParamGetShortVector(lowerModel)

    upperModel <- do.call(PCMParamUpperLimit, c(list(model), argsPCMParamUpperLimit))
    upperVecParams <- PCMParamGetShortVector(upperModel)

    config <- do.call(configOptimAndMCMC,
                      c(list(lik = lik,
                             parLower = lowerVecParams,
                             parUpper = upperVecParams),
                        argsConfigOptimAndMCMC) )
  }

  res <- runOptimAndMCMC(lik, prior, input.data, config, verbose)

  res$X <- X
  res$tree <- tree
  res$modelInit <- model
  res$SE <- SE
  res$lik <- lik
  res$config <- config
  res$lowerVecParams <- lowerVecParams
  res$lowerModel <- lowerModel
  res$upperVecParams <- upperVecParams
  res$upperModel <- upperModel
  res$PCMOptions <- PCMOptions()

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

  attr(value, "df") <- PCMParamCount(object$modelInit,
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

  model <- object$modelInit

  if(!is.null(object$modelOptim)) {
    model <- object$modelOptim
  }

  par <- PCMGetVecParamsRegimesAndModels(model, object$tree, ...)
  c(par, numParam = PCMParamCount(model))
}

