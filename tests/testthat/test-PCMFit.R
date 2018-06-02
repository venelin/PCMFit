source("GenerateTestTreeAndData.R")
fit2 <- PCMFit(values, tree.ab.singles, model, metaI = metaInfo,
               argsConfigOptimAndMCMC = list(nCallsOptim = 20, genInitVerbose = TRUE),
               verbose = TRUE)

test_that("ML value is better than original",
          expect_true(logLik(fit2) < 0 && logLik(fit2) > likOrig))

