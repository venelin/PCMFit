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
#' @importFrom PCMBase PCMCreateLikelihood PCMInfo PCMParamCount PCMParamGetShortVector PCMParamLoadOrStore PCMParamLowerLimit PCMParamUpperLimit PCMParamRandomVecParams PCMOptions
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
  verbose = FALSE) {

  if(is.function(metaI)) {
    metaI <- metaI(X = X, tree = tree, model = model, SE = SE)
  }

  lik <- PCMCreateLikelihood(
    X = X, tree = tree, model = model, SE = SE, metaI = metaI,
    positiveValueGuard = positiveValueGuard)

  lowerModel <- do.call(PCMParamLowerLimit, c(list(model), argsPCMParamLowerLimit))
  lowerVecParams <- PCMParamGetShortVector(lowerModel)

  upperModel <- do.call(PCMParamUpperLimit, c(list(model), argsPCMParamUpperLimit))
  upperVecParams <- PCMParamGetShortVector(upperModel)


  matParInitRunif <- PCMParamRandomVecParams(
    o = model,
    k = PCMNumTraits(model),
    R = PCMNumRegimes(model),
    n = numRunifInitVecParams,
    argsPCMParamLowerLimit = argsPCMParamLowerLimit,
    argsPCMParamUpperLimit = argsPCMParamUpperLimit
  )
  matParInitGuess <- GuessInitVecParams(
    o = model,
    k = PCMNumTraits(model),
    R = PCMNumRegimes(model),
    n = numGuessInitVecParams,
    argsPCMParamLowerLimit = argsPCMParamLowerLimit,
    argsPCMParamUpperLimit = argsPCMParamUpperLimit,
    X = X, tree = tree, SE = SE, varyTheta = FALSE)
  matParInitGuessVaryTheta <- GuessInitVecParams(
    o = model,
    k = PCMNumTraits(model),
    R = PCMNumRegimes(model),
    n = numGuessInitVecParams,
    argsPCMParamLowerLimit = argsPCMParamLowerLimit,
    argsPCMParamUpperLimit = argsPCMParamUpperLimit,
    X = X, tree = tree, SE = SE, varyTheta = TRUE)

  matParInit <- rbind(
    matParInit,
    matParInitRunif,
    matParInitGuess,
    matParInitGuessVaryTheta)

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
      cat("Evaluating likelihood at ", nrow(matParInit), " parameter vectors...")
    }
    valParInitOptim <- apply(matParInit, 1, lik)

    if(verbose) {
      cat(
        "Taking the top-", numCallsOptim,
        " parameter combinations sorted by decreasing log-likelihood value...")
    }

    topVal <- order(valParInitOptim, decreasing = TRUE)[1:numCallsOptim]
    matParInit <- matParInit[topVal,, drop=FALSE]
  }

  res <- as.list(environment())
  # These objects tend to be very big. the lik function and the metaI object
  # can be recreated.
  res$lik <- res$metaI <- res$matParInit <-
    res$matParInitRunif <- res$matParInitGuess <-
    res$matParInitGuessVaryTheta <- NULL

  res$PCMOptions <- PCMOptions()

  res <- c(
    res,
    runOptim(
      lik = lik,
      parLower = lowerVecParams,
      parUpper = upperVecParams,
      matParInit = matParInit,
      control = control,
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

#'@export
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

# A utility function used to save the maximum point of f (used as wrapper for
# likelihod-functions).
memoiseMax <- function(f, par, memo, verbose) {
  countMemo <- mget('count', envir = memo, ifnotfound = list(0))$count
  valMemo <- mget('val', envir = memo, ifnotfound = list(-Inf))$val
  valPrev <- mget('valPrev', envir = memo, ifnotfound = list(-Inf))$valPrev
  valDelta <- mget('valDelta', envir = memo, ifnotfound = list(NA))$valDelta
  parMemo <- mget('par', envir = memo, ifnotfound = list(NULL))$par

  assign("count", countMemo + 1, pos = memo)

  val <- f(par)

  # if the returned val higher than the current max-value, store this val:
  if(valMemo < val) {
    assign('par', par, pos = memo)
    assign('val', val, pos = memo)
    assign('valDelta', val - valMemo, pos = memo)

    if(verbose) {
      cat('\nCall ', countMemo, ': value on par=(',
          toString(round(par, 6)), "): ", val, "\n", sep = "")
    }
  }

  # store the returned val in memo
  assign('valPrev', val, pos = memo)

  val
}


#' Calling optim a number of times
#' @param lik a function of numeric vector argument. Possible signatures are
#' function(par) {} or function(par, input) {}. This argument is mandatory and
#' doesn't have a  default value.
#' @param parLower,parUpper numeric vectors of equal length p: the number of
#' parameters of the function lik;
#' @param matParInit a numeric matrix of p columns.
#' @param control a list passed to optim().
#' @param verbose logical indicating whether informative messages should be
#' printed on the console during the run.
#'
#' @return a named list
#' @seealso \code{\link{configOptim}}.
#'
#' @importFrom stats optim
#' @importFrom foreach foreach %do% %dopar%
#' @import data.table
runOptim <- function(
  lik,
  parLower,
  parUpper,
  matParInit,
  control = NULL,
  verbose = TRUE) {

  res <- list(Optim = NULL)
  class(res) <- "ResultOptim"

  tryCatch({

    if(!(is.vector(parLower) && is.vector(parUpper) &&
         length(parLower) == length(parUpper) &&
         length(parLower) > 0 &&
         ncol(matParInit) == length(parLower))) {
      stop(
        paste(
          "parLower, parUpper should be numeric vectors of non-zero length, p;",
          " matParInit should be a numeric columns with p columns."))
    }

    if(!isTRUE(all(parLower < parUpper))) {
      stop("All elements of parLower should be smaller than the corresponding elements in parUpper.")
    }

    if(!is.function(lik)) {
      stop("Expecting lik to be a function(par).")
    }

    memoMaxLoglik <- new.env()
    fnForOptim <- function(par) {
      memoiseMax(lik, par = par, memoMaxLoglik, verbose)
    }

    listCallsOptim <- list()

    for(iOptimTry in seq_len(nrow(matParInit))) {

      listCallsOptim[[iOptimTry]] <- list()

      listCallsOptim[[iOptimTry]]$parStart <- parInit <- matParInit[iOptimTry, ]
      listCallsOptim[[iOptimTry]]$valueStart <- fnForOptim(parInit)

      if(!(isTRUE(all(parInit >= parLower) &&
                  all(parInit <= parUpper)))) {
        # this should in principle never happen, because parInit is not user-specified.
        #
        warning(
          paste("Skipping optim try #", iOptimTry, ":",
                "All parameters in parInit should be between \n parLower=c(",
                toString(parLower), ") # and \n parUpper=c(",
                toString(parUpper), ") #, but were \n parInit = c(",
                toString(parInit), ")"))
        next
      }

      if(is.null(control)) {
        control <- list()
      } else {
        control <- as.list(control)
      }

      # ensure that optim does maximization.
      control$fnscale <- -1

      if(length(parInit) > 0) {
        res.optim.call <-
          optim(fn = fnForOptim,
                par = parInit, lower = parLower, upper = parUpper,
                method = 'L-BFGS-B', control = control)

        listCallsOptim[[iOptimTry]]$parEnd <- res.optim.call$par
        listCallsOptim[[iOptimTry]]$valueEnd <- res.optim.call$value
        listCallsOptim[[iOptimTry]]$counts <- res.optim.call$counts
        listCallsOptim[[iOptimTry]]$convergence <- res.optim.call$convergence

        if(verbose) {
          cat("\nCall to optim no.", iOptimTry,
              ": starting from ",
              toString(round(parInit, 4)), ": ",
              round(listCallsOptim[[iOptimTry]]$valueStart, 4), "\n",
              "parLower = c(", toString(round(parLower, 4)), ")\n",
              "parUpper = c(", toString(round(parUpper, 4)), ")\n",
              "value: ", round(res.optim.call$value, 4), "\n",
              "par: ", toString(round(res.optim.call$par, 4)), "\n",
              "convergence: ", res.optim.call$convergence, "\n",
              "counts: ", toString(res.optim.call$counts), "\n",
              "message: ", toString(res.optim.call$message), "\n")
        }
      } else {
        stop("optim: parameter vector has zero length.")
      }
    }

    maxPar <- get("par", pos = memoMaxLoglik)
    maxValue <- get("val", pos = memoMaxLoglik)
    callCount <- get("count", pos = memoMaxLoglik)

    res.optim <- list(par = maxPar, value = maxValue, count = callCount,
                      listCallsOptim = listCallsOptim)

    res$Optim <- res.optim

  },
  interrupt = function() {
    return(res)
  },
  error = function(e) {
    cat("Error in runOptim:", toString(e), "trace:")
    traceback()
    stop(e)
  })

  res
}
