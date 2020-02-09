# A utility function used to save the maximum point of f (used as wrapper for
# likelihod-functions).
MemoiseMax <- function(f, par, memo, verbose) {
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
#' @inheritParams PCMFit
#' @param parLower,parUpper numeric vectors of equal length p: the number of
#' parameters of the function lik;
#' @param matParInit a numeric matrix of p columns.
#' @param control a list passed to optim().
#' @param doParallel logical indicating if optim calls should be executed in
#' parallel using the \code{foreach() \%dopar\% {}} construct. Default: FALSE.
#' @param verbose logical indicating whether informative messages should be
#' printed on the console during the run.
#'
#' @return a named list
#'
#' @importFrom stats optim
#' @importFrom foreach foreach %do% %dopar%
#' @import data.table
RunOptim <- function(
  X, tree, model, SE, metaI, positiveValueGuard,
  parLower,
  parUpper,
  matParInit,
  control = NULL,
  doParallel = FALSE,
  verbose = TRUE) {

  # prevent 'no visible binding' notes
  iOptimTry <- NULL

  res <- list(Optim = NULL)
  class(res) <- "ResultOptim"

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

  `%op%` <- if(isTRUE(doParallel) ||
               (is.numeric(doParallel) && doParallel > 1)) {
    `%dopar%`
  } else {
    `%do%`
  }

  globalFuns <- unname(unlist(
    sapply(ls(.GlobalEnv), function(n) if(is.function(.GlobalEnv[[n]])) n else NULL)))


  listCallsOptim <- foreach(iOptimTry = seq_len(nrow(matParInit)),
                            .export = globalFuns) %op% {
                              memoMaxLoglik <- new.env()

                              lik <- PCMCreateLikelihood(
                                X = X, tree = tree, model = model, SE = SE,
                                metaI = metaI,
                                positiveValueGuard = positiveValueGuard)

                              fnForOptim <- function(par) {
                                MemoiseMax(lik, par = par, memoMaxLoglik, verbose)
                              }

                              resIter <- list()

                              resIter$parStart <- parInit <- matParInit[iOptimTry, ]
                              resIter$valueStart <- fnForOptim(parInit)

                              if(!(isTRUE(all(parInit >= parLower) &&
                                          all(parInit <= parUpper)))) {
                                # this should in principle never happen, because parInit is not user-specified.
                                #
                                if(verbose) {
                                  cat(
                                    paste("Skipping optim try #", iOptimTry, ":",
                                          "All parameters in parInit should be between \n parLower=c(",
                                          toString(parLower), ") # and \n parUpper=c(",
                                          toString(parUpper), ") #, but were \n parInit = c(",
                                          toString(parInit), ")"))
                                }
                                resIter$maxPar <- parInit
                                resIter$maxValue <- getOption("PCMBase.Value.NA", default = -1e+20)
                                resIter$callCount <- 0L

                              } else {
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

                                  resIter$parEnd <- res.optim.call$par
                                  resIter$valueEnd <- res.optim.call$value
                                  resIter$counts <- res.optim.call$counts
                                  resIter$convergence <- res.optim.call$convergence

                                  if(verbose) {
                                    cat("\nCall to optim no.", iOptimTry,
                                        ": starting from ",
                                        toString(round(parInit, 4)), ": ",
                                        round(resIter$valueStart, 4), "\n",
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

                                resIter$maxPar <- get("par", pos = memoMaxLoglik)
                                resIter$maxValue <- get("val", pos = memoMaxLoglik)
                                resIter$callCount <- get("count", pos = memoMaxLoglik)
                              }

                              resIter
                            }

  iBestCallOptim <- which.max(sapply(listCallsOptim, function(.) .$maxValue))
  res.optim <- list(
    par = listCallsOptim[[iBestCallOptim]]$maxPar,
    value = listCallsOptim[[iBestCallOptim]]$maxValue,
    count = listCallsOptim[[iBestCallOptim]]$callCount,
    listCallsOptim = listCallsOptim)

  res$Optim <- res.optim

  res
}
