# Copyright 2016-2020 Venelin Mitov
#
# This file is part of PCMFit.
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

#' A \code{\link{data.table}} of \code{\link{MGPM}} objects
#'
#' @param o an object, such as a character string denoting a path to a '.RData',
#' '.rda' file where an object named as indicated by \code{objectName} is stored,
#' or a \code{\link{MGPMEnvironment}} object, or a \code{\link{MGPM}}, or a list of
#' \code{MGPM} objects or a list of \code{double} vectors (see argument \code{env}).
#' @param objectName \code{NULL} or a character string denoting an R object stored
#' in the file indicated by \code{o}. This argument is ignored if \code{o} is not a
#' character string (file-path).
#' @param env a \code{\link{MGPMEnvironment}} object. This must be specified only if
#' \code{o} is a list of \code{double} vectors. Otherwise this argument is ignored.
#' @param treeEDExpr a character string denoting a ED-expression
#' (see \code{\link{PCMBase::PCMTreeEvalNestedEDxOnTree}}). Like the \code{env}
#' argument, this must be specified only if \code{o} is a list of \code{double}
#' vectors. Otherwise this argument is ignored.
#' @param verbose logical (default \code{FALSE}) indicating if some log-info
#' should be printed to the console.
#'
#' @details This is a S3 generic function. See \code{methods(MGPMTable)} for
#' implementing methods.
#' @return a \code{MGPMTable} object.
#' @export
MGPMTable <- function(
  o, objectName = NULL, env = NULL, treeEDExpr = "tree", verbose = FALSE) {

  UseMethod("MGPMTable", o)
}

#' @export
MGPMTable.default <- function(
  o, objectName = NULL, env = NULL, treeEDExpr = "tree", verbose = FALSE) {

  stop(
    paste0("MGPMTable.default:: Don't know how to convert o to a MGPMTable; class(o)=",
           toString(o), "."))
}

#' @export
MGPMTable.MGPMTable <- function(
  o, objectName = NULL, env = NULL, treeEDExpr = "tree", verbose = FALSE) {
  o
}

#' @export
MGPMTable.character <- function(
  o, objectName = NULL, env = NULL, treeEDExpr = "tree", verbose = FALSE) {

  if(!file.exists(o)) {
    stop(paste0("MGPMTable.character:: file ", o, " does not exist."))
  } else if(!is.character(objectName)) {
    stop(paste0("MGPMTable.character:: objectName should be a character denoting the ",
                "name of an object stored in the file ", o, "."))
  } else {
    load(o, verbose = verbose)
    object <- try(get(objectName), silent = TRUE)
    if(inherits(object, "try-error")) {
      stop(paste0("MGPMTable.character:: could not locate an object ", objectName,
                  " in file ", o, "."))
    } else {
      if(verbose) {
        cat("Calling MGPMTable on loaded object...")
        MGPMTable(object, objectName = NULL, env = env, treeEDExpr = treeEDExpr,
                  verbose = verbose)
      }
    }
  }
}

#' @export
MGPMTable.MGPMEnvironment <- function(
  o, objectName = NULL, env = NULL, treeEDExpr = "tree", verbose = FALSE) {

  dt <- data.table(
    hash.tree = character(),
    hash.n = character(),
    hash.l = character(),
    hash.r = character(),
    hash.m = character(),
    hash.v = character(),

    mgpm = list(),

    ll = double(),
    df = integer(),
    nobs = integer(),
    score = double())
  class(dt) <- c("MGPMTable", class(dt))
  setkey(dt, hash.tree, hash.n, hash.l, hash.r, hash.m, hash.v)
  dt
}

#' @export
MGPMTable.MGPM <- function(
  o, objectName = NULL, env = NULL, treeEDExpr = "tree", verbose = FALSE) {

  dt <- data.table(
    hash.tree = character(),
    hash.n = character(),
    hash.l = character(),
    hash.r = character(),
    hash.m = character(),
    hash.v = character(),

    mgpm = list(),

    ll = double(),
    df = integer(),
    nobs = integer(),
    score = double())
  class(dt) <- c("MGPMTable", class(dt))
  setkey(dt, hash.tree, hash.n, hash.l, hash.r, hash.m, hash.v)
  dt
}

#' Check if an object is a MGPMTable.
#' @param o an R object.
#' @return a logical
#' @export
is.MGPMTable <- function(o) {
  inherits(o, "MGPMTable")
}

#   # prevent 'no visible binding' notes
#   hashCodeTree <- hashCodeStartingNodesRegimesLabels <-
#     hashCodeMapping <- NULL
#
#   tableFits <- tableFitsPrev
#
#   if(!is.null(fitMappingsPrev)) {
#     tableFits <- tableFitsPrev <- fitMappingsPrev$tableFits
#     if(!identical(modelTypes, fitMappingsPrev$arguments$modelTypes)) {
#
#       # this should remap the model-type indices in the fit vectors, show
#       # table fits is correctly converted.
#       tableFits <- RetrieveFittedModelsFromFitVectors(
#         fitMappings = fitMappingsPrev, tableFits = tableFitsPrev,
#         modelTypesNew = modelTypes)
#     }
#   }
#
#   if(!is.data.table(tableFits)) {
#     if(verbose) {
#       cat("Initiating tableFits...\n")
#     }
#     tableFits = data.table(hashCodeTree = character(),
#                            hashCodeStartingNodesRegimesLabels = character(),
#                            hashCodeMapping = character(),
#                            treeEDExpression = character(),
#                            startingNodesRegimesLabels = list(),
#                            mapping = list(),
#                            fitVector = list(),
#                            logLik = double(),
#                            df = integer(),
#                            nobs = integer(),
#                            score = double(),
#                            duplicated = logical())
#   }
#
#   setkey(tableFits, hashCodeTree,hashCodeStartingNodesRegimesLabels,hashCodeMapping )
#   attr(tableFits, "modelTypes") <- modelTypes
#   tableFits
# }
