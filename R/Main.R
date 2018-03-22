# Copyright 2018 Venelin Mitov
#
# This file is part of PCMStep
#
# PCMStep is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PCMStep is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PCMBase.  If not, see <http://www.gnu.org/licenses/>.


#' Generate the next mapping of a set of models to a number of regimes
#' @param mapping a character vector with elements from models
#' (repetitions are allowed) denoting a current model-to-regime mapping.
#' @param models a character vector with the class-names of PCM models
#' @return a character vector of the same length as mapping with elements from models.
#' @export
PCMNextMapping <- function(mapping, models) {
  R <- length(mapping)
  mappingInd <- match(mapping, models)
  if(any(is.na(mappingInd))) {
    print(mapping)
    print(models)
    print(mappingInd)
    stop(paste0("ERR:04000:PCMStep:Main.R:PCMStep:: mapping should have
         length ", R, " and contain only elements among models"))
  }


  numModels <- length(models)

  carry <- 1
  for(pos in R:1) {
    if(mappingInd[pos] + carry <= numModels) {
      mappingInd[pos] <- mappingInd[pos] + carry
      carry <- 0
    } else {
      mappingInd[pos] <- 1
    }
  }
  res <- models[mappingInd]
  names(res) <- names(mapping)
  res
}

#' Iterator over combinations with repetions of a given set of models
#' @param mapping a vector of elements from models giving the initial combination
#' @param models a vector of unique elements to choose from when building the
#' combinations.
#' @return an iterator object with S3 class c("imapping", "abstractiter", "iter").
#' Calling repeatedly nextElem on this object iterates over all possible combinations
#' with repetitions.
#' @examples
#' it <- imapping(c("BM3", "BM3"), c("BM3", "OU3", "JOU3"))
#' nextElem(it)
#' nextElem(it)
#' @export
imapping <- function(mapping, models) {
  state <- new.env()

  state$initial <- mapping
  state$current <- NULL # at initial state before calling nextEl

  nextEl <- function() {
    if(is.null(state$current)) {
      state$current <- state$initial
      state$current
    } else {
      nextMapping <- NextMapping(state$current, models)
      if(isTRUE(all(state$initial==nextMapping))) {
        stop("StopIteration", call. = FALSE)
      } else {
        state$current <- nextMapping
        nextMapping
      }
    }
  }
  obj <- list(nextElem=nextEl)
  class(obj) <- c('imapping', 'abstractiter', 'iter')
  obj
}

#' Plot a tree with regimes
#' @importFrom data.table data.table
#' @importFrom grDevices hcl
#' @importFrom PCMBase PCMNumUniqueRegimesTree PCMUniqueRegimesTree
#' @importFrom ggtree ggtree %<+%
#' @importFrom ggplot2 aes scale_color_manual
#' @export
PCMPlotTree <- function(tree) {
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  N <- PCMNumTips(tree)
  R <- PCMNumUniqueRegimesTree(tree)
  palette <- gg_color_hue(R)
  names(palette) <- PCMUniqueRegimesTree(tree)

  data <- rbind(data.table(node = tree$edge[, 2], regime = tree$edge.regime),
                data.table(node = N+1, regime = NA))

  plotTree <- ggtree(tree, layout = 'fan', open.angle = 8, size=.25) %<+% data

  plotTree + aes(color = regime) +
    scale_color_manual(name = "regime", values = palette)
}
