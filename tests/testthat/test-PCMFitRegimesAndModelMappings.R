source("GenerateTestTreeAndData.R")

modelTypes <-  c("BM__noX0__noSigmae_x", "OU__noX0__noSigmae_x")
argsMRG <- list(X0 = list(default = rep(as.double(0), 2), type = c("gvector", "full")),
                Sigmae_x = list(default = matrix(0, 2, 2), type = c("gmatrix", "fixed"),
                                description = "Fixed upper triangular Choleski factor of the variance-covariance matrix for the non-phylogenetic trait component"))

argsPCMLowerBound <- list(lowerBoundValue = -10, lowerBoundValuePositiveDiag = 0)
argsPCMUpperBound <- list(upperBoundValue = 10, upperBoundValuePositiveDiag = 10)

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)

tableFits <- NULL

# tableFits <- rbindlist(lapply(list.files(".", pattern="fits.*.RData"), function(f)  {
#   load(f)
#   data
# }))
#
# print(tableFits)
# tableFits[, duplicated:=FALSE]
#
# if(nrow(tableFits) == 0) {
#   tableFits <- NULL
# }

if(!exists("cluster") || is.null(cluster)) {
  cluster <- parallel::makeCluster(parallel::detectCores(logical = TRUE), outfile = "log.txt")
  doParallel::registerDoParallel(cluster)
}

tableFits <- PCMFitRegimesAndModelMappingsRecursive(
  values, tree.ab.singles, modelTypes = modelTypes, metaIFun = PCMInfoCpp, positiveValueGuard = 1000,


  tableFits = tableFits,

  minCladeSize = 35,

  argsMRG = argsMRG,
  argsPCMLowerBound = argsPCMLowerBound,
  argsPCMUpperBound = argsPCMUpperBound,
  argsConfigOptimAndMCMC1 = list(nCallsOptim = 10, genInitNumEvals = 100, genInitVerbose = FALSE),
  argsConfigOptimAndMCMC2 = list(nCallsOptim = 10, genInitNumEvals = 10, genInitVerbose = FALSE),

  numJitterAllRegimeFits = 100, numJitterRootRegimeFit = 100,

  printFitVectorsToConsole = TRUE,
  doParallel = TRUE,
  verbose = TRUE)


# Don't forget to destroy the parallel cluster to avoid leaving zombie worker-processes.
parallel::stopCluster(cluster)
cluster <- NULL




