library(testthat)
context("PCMParam")

library(PCMBase)
library(PCMFit)

load("testobjects.RData")

set.seed(1)

k <- PCMNumTraits(model.ab.123)
R <- PCMNumRegimes(model.ab.123)

randVecs1 <- PCMParamRandomVecParams(o = model.ab.123, k = k, R = R, n = 10)

randVecs3 <- guessInitVecParams(
  o = model.ab.123, k = k, R = R, n = 10,
  X = traits.ab.123[, seq_len(PCMTreeNumTips(tree.ab))],
  tree = tree.ab)



