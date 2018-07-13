source("GenerateTestTreeAndData.R")

library(PCMFit)
if(require(PCMBaseCpp)) {
  metaInfo1 <- PCMInfoCpp(values, tree.a, model1)
  metaInfo2 <- PCMInfoCpp(values, tree.a, model2)
} else {
  metaInfo1 <- PCMInfo(values, tree.a, model1)
  metaInfo2 <- PCMInfo(values, tree.a, model2)
}

fit1 <- PCMFit(values, tree.a, model1, metaI = metaInfo1,
                 argsConfigOptimAndMCMC = list(genInitNumEvals = 500000, nCallsOptim = 500, genInitVerbose = TRUE),
                 verbose = TRUE)


fit2 <- PCMFit(values, tree.a, model2, metaI = metaInfo2,
               argsConfigOptimAndMCMC = list(genInitNumEvals = 100000, nCallsOptim = 100, genInitVerbose = TRUE),
               verbose = TRUE)


test_that("ML value is better than original",
          expect_true(logLik(fit1) < 0 && logLik(fit1) > likOrig))


test_that("fit1 outperforms fit2",
          expect_true(AIC(fit2) > AIC(fit1)))

