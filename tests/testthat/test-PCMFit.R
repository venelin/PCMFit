library(testthat)

context("PCMFit")

library(PCMFit)
library(PCMBase)
library(PCMBaseCpp)

load("testobjects.RData")

options(PCMBase.Value.NA = -1e20)

PCMLik(X = traits.a.1, tree = tree.a, model = model.a.1)

fit.a.1 <- PCMFit(
  X = traits.a.1,
  tree = tree.a,
  model = model.a.1,
  metaI = PCMInfoCpp,
  verbose = TRUE)

fit.a.1.fromTrue <- PCMFit(
  X = traits.a.1,
  tree = tree.a,
  model = model.a.1,
  metaI = PCMInfoCpp,
  matParInit = matrix(PCMParamGetShortVector(model.a.1,
                                             k = PCMNumTraits(model.a.1),
                                             R = PCMNumRegimes(model.a.1)), 1),
  verbose = TRUE)




PCMLik(X = traits.a.123, tree = tree.a, model = model.a.123)

fit.a.123 <- PCMFit(
  X = traits.a.123,
  tree = tree.a,
  model = model.a.123,
  metaI = PCMInfoCpp,
  verbose = TRUE)

fit.a.123.fromTrue <- PCMFit(
  X = traits.a.123,
  tree = tree.a,
  model = model.a.123,
  metaI = PCMInfoCpp,
  matParInit = matrix(PCMParamGetShortVector(model.a.123,
                                             k = PCMNumTraits(model.a.123),
                                             R = PCMNumRegimes(model.a.123)), 1),
  verbose = TRUE)


logLik(fit.a.1)
