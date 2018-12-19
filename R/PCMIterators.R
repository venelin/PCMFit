#' Generate the next mapping of a set of model-types to a number of regimes
#' @param mapping a vector with elements from \code{models}.
#' (repetitions are allowed) denoting a current model-to-regime mapping.
#' @param modelTypes a vector of PCM model-types or integer indices.
#' @return a vector of the same length as mapping with elements from \code{modelTypes}.
#' @export
PCMNextMapping <- function(mapping, modelTypes) {
  R <- length(mapping)
  mappingInd <- match(mapping, modelTypes)
  if(any(is.na(mappingInd))) {
    stop(paste0("ERR:04100:PCMFit:PCMFit.R:PCMNextMapping:: mapping should have
                length ", R,
                " and contain only elements among modelTypes; mapping = (",
                toString(mapping),")",
                ", modelTypes=(", toString(modelTypes), "), mappingInd=(",
                toString(mappingInd), ")."))
  }


  numModels <- length(modelTypes)

  carry <- 1
  for(pos in R:1) {
    if(mappingInd[pos] + carry <= numModels) {
      mappingInd[pos] <- mappingInd[pos] + carry
      carry <- 0
    } else {
      mappingInd[pos] <- 1
    }
  }
  res <- modelTypes[mappingInd]
  names(res) <- names(mapping)
  res
}

#' Iterator over combinations with repetions of a given set of modelTypes
#' @param mapping a vector of elements from modelTypes giving the initial combination
#' @param modelTypes a vector of unique elements to choose from when building the
#' combinations.
#' @return an iterator object with S3 class c("imapping", "abstractiter", "iter").
#' Calling repeatedly nextElem on this object iterates over all possible combinations
#' with repetitions of the same length as the argument \code{mapping}.
#' @examples
#' it <- PCMIteratorMapping(c("BM", "BM"), c("BM", "OU", "JOU"))
#' iterators::nextElem(it)
#' iterators::nextElem(it)
#' @export
PCMIteratorMapping <- function(mapping, modelTypes) {
  state <- new.env()

  state$initial <- mapping
  state$current <- NULL # at initial state before calling nextEl

  nextEl <- function() {
    if(is.null(state$current)) {
      state$current <- state$initial
      state$current
    } else {
      nextMapping <- PCMNextMapping(state$current, modelTypes)
      if(isTRUE(all(state$initial==nextMapping))) {
        stop("StopIteration", call. = FALSE)
      } else {
        state$current <- nextMapping
        nextMapping
      }
    }
  }
  obj <- list(nextElem=nextEl)
  class(obj) <- c('PCMIteratorMapping', 'abstractiter', 'iter')
  obj
}


#' Generate the next mapping of model-types chosen from regime-specific sets of possible model-types.
#' @param mapping a vector with elements from \code{modelTypes}.
#' (repetitions are allowed) denoting a current model-to-regime mapping.
#' @param modelTypes a vector of unique model-types, e.g. \code{c("BM", "OU")}.
#' @param allowedModelTypesIndices a list of the same length as \code{mapping} with integer vector elements
#' or NULLs. When an element of this list is an integer vector its elements denote unique positions in
#' modelTypes, i.e. the allowed model-types for the regime at that position in mapping.
#' @return a vector of the same length as mapping with elements from \code{modelTypes}.
#' @export
PCMNextMapping2 <- function(mapping, modelTypes, allowedModelTypesIndices) {
  R <- length(mapping)
  numModels <- length(modelTypes)

  if(length(allowedModelTypesIndices) != length(mapping)) {
    stop(paste0("ERR:04130:PCMFit:PCMFit.R:PCMNextMapping2:: mapping and allowedModelTypesIndices should
                be the same length, ", R, "."))
  }

  mappingInd <- match(mapping, modelTypes)
  #cat("mappingInd: ", toString(mappingInd), "\n")
  mappingInd2 <- sapply(1:R, function(pos) {
    if(!is.null(allowedModelTypesIndices[[pos]])) {
      #cat("pos = ", pos, "; allowedModelTypesIndicies[[pos]] = ", toString(allowedModelTypesIndices[[pos]]), "\n")
      #cat("pos = ", pos, "; mappingInd[pos] = ", mappingInd[pos], "\n")
      match(mappingInd[pos], allowedModelTypesIndices[[pos]])
    } else {
      # allowedModelTypesIndices[[pos]] == NULL means that all model types are allowed at this position
      match(mappingInd[pos], 1:numModels)
    }
  })


  if(any(is.na(mappingInd2))) {
    stop(paste0("ERR:04100:PCMFit:PCMFit.R:PCMNextMapping2:: mapping should have
                length ", R, " allowedModelTypesIndices; mapping = (", toString(mapping),")",
                ", modelTypes=(", toString(modelTypes), "), mappingInd2=(", toString(mappingInd2), ")."))
  }


  numsAllowedModelTypes <- sapply(1:R, function(pos) {
    if(is.null(allowedModelTypesIndices[[pos]])) {
      numModels
    } else {
      length(allowedModelTypesIndices[[pos]])
    }
  })


  carry <- 1
  for(pos in R:1) {
    if(mappingInd2[pos] + carry <= numsAllowedModelTypes[pos]) {
      mappingInd2[pos] <- mappingInd2[pos] + carry
      carry <- 0
    } else {
      mappingInd2[pos] <- 1
    }
  }


  #cat("mappingInd2=", toString(mappingInd2), "\n")
  res <- sapply(1:R, function(pos) {
    if(is.null(allowedModelTypesIndices[[pos]])) {
      modelTypes[mappingInd2[pos]]
    } else {
      modelTypes[allowedModelTypesIndices[[pos]][mappingInd2[pos]]]
    }
  })

  #cat("indices=", toString(indices), "\n")
  #res <- modelTypes[indices]
  names(res) <- names(mapping)
  res
  }

#' Iterator over combinations with repetions of a given set of modelTypes
#' @param mapping a vector of elements from modelTypes giving the initial combination
#' @param modelTypes a vector of unique elements to choose from when building the
#' combinations.
#' @param allowedModelTypesIndices a list of the same length as \code{mapping} with integer vector elements
#' or NAs. When an element of this list is an integer vector its elements denote unique positions in
#' modelTypes, i.e. the allowed model-types for the regime at that position in mapping.
#' @return an iterator object with S3 class c("imapping", "abstractiter", "iter").
#' Calling repeatedly nextElem on this object iterates over all possible combinations
#' with repetitions of the same length as the argument \code{mapping}.
#' @examples
#' it <- PCMIteratorMapping(c("BM", "BM"), c("BM", "OU", "JOU"))
#' iterators::nextElem(it)
#' iterators::nextElem(it)
#' @export
PCMIteratorMapping2 <- function(mapping, modelTypes, allowedModelTypesIndices) {
  state <- new.env()

  state$initial <- mapping
  state$current <- NULL # at initial state before calling nextEl

  nextEl <- function() {
    if(is.null(state$current)) {
      state$current <- state$initial
      state$current
    } else {
      nextMapping <- PCMNextMapping2(state$current, modelTypes, allowedModelTypesIndices)
      if(isTRUE(all(state$initial==nextMapping))) {
        stop("StopIteration", call. = FALSE)
      } else {
        state$current <- nextMapping
        nextMapping
      }
    }
  }
  obj <- list(nextElem=nextEl)
  class(obj) <- c('PCMIteratorMapping2', 'abstractiter', 'iter')
  obj
}
