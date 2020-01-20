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

#' Environment of a MGPM
#'
#' A \code{\link{MGPMEnvironment}} represents an R environment containing the
#' necessary data and meta-information that remains constant (e.g. candidate model
#' types for different MGPM regimes) or change rarely (e.g. metaI objects) during a
#' MGPM inference.
#'
#' @inheritParams PCMBase::PCMLik
#' @param model a MixedGaussian model used as a template to build MixedGaussian
#' model objects.
#' @param tipModelType NULL or a character string indicating a member of
#' \code{model}. This argument allows to set a specific regime for all terminal
#' branches, such as a white noise model.
#' @param parent a parent environment. Default: \code{parent.frame()}.
#' @param ... any additional named objects to be included in the MGPMEnvironment
#' environment.
#'
#' @return an object of S3 class 'MGPMEnvironment'.
#' @seealso \code{\link{MGPM}}
#' @export
MGPMEnvironment <- function(
  X, tree, model, SE = matrix(0, PCMNumTraits(model), PCMTreeNumTips(tree)),
  tipModelType = NULL, parent = parent.frame(), ...) {

  env <- list2env(list(
    k = nrow(X),
    X = X,
    SE = SE,
    tipModelType = tipModelType,
    modelTemplate = model,
    treeOriginal = tree,
    tree = PCMTree(tree),
    options = PCMOptions(),
    ...),
    parent = parent)

  PCMTreeSetLabels(env$tree)

  class(env) <- "MGPMEnvironment"

  env
}

#' Check if an object is of S3 class 'MGPMEnvironment'
#' @param o an object.
#' @return a logical.
#' @export
is.MGPMEnvironment <- function(o) {
  inherits(o, "MGPMEnvironment")
}
