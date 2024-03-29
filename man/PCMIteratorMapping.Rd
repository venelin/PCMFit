% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCMIterators.R
\name{PCMIteratorMapping}
\alias{PCMIteratorMapping}
\title{Iterator over combinations with repetions of a given set of modelTypes}
\usage{
PCMIteratorMapping(
  mapping,
  modelTypes,
  allowedModelTypesIndices = rep(list(NULL), length(mapping))
)
}
\arguments{
\item{mapping}{a vector of elements from modelTypes giving the initial combination}

\item{modelTypes}{a vector of unique elements to choose from when building the
combinations.}

\item{allowedModelTypesIndices}{a list of the same length as \code{mapping}
with integer vector elements or NULLs. When an element of this list is an
integer vector its elements denote unique positions in modelTypes, i.e. the
allowed model-types for the regime at that position in mapping. Default
value: \code{rep(list(NULL), length(mapping))}.}
}
\value{
an iterator object with S3 class \code{c("imapping", "abstractiter",
"iter")}.
Calling repeatedly nextElem on this object iterates over all possible
combinations with repetitions of the same length as the argument
\code{mapping}.
}
\description{
Iterator over combinations with repetions of a given set of modelTypes
}
\examples{
library(iterators)
it <- PCMIteratorMapping(c(1, 1), c(1, 2, 3))
nextElem(it)
nextElem(it)
nextElem(it)
nextElem(it)

it <- PCMIteratorMapping(c(1, 1), c(1, 2, 3), list(NULL, 1:2))
nextElem(it)
nextElem(it)
nextElem(it)
nextElem(it)
nextElem(it)
nextElem(it)

it <- PCMIteratorMapping(c("BM", "BM"), c("BM", "OU", "JOU"))
nextElem(it)
nextElem(it)
}
