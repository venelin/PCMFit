library(data.table)
library(PCMFit)
prefixFiles <- "local_201812171_"

if(!exists("cluster") || is.null(cluster)) {
  if(require(doMPI)) {
    # using MPI cluster as distributed node cluster (possibly running on a cluster)
    # Get the number of cores. Assume this is run in a batch job.
    p = strtoi(Sys.getenv('LSB_DJOB_NUMPROC'))
    cluster <- startMPIcluster(count = p-1, verbose = TRUE)
    doMPI::registerDoMPI(cluster)
  } else {
    cluster <- parallel::makeCluster(parallel::detectCores(logical = TRUE),
                                     outfile = paste0("log_", prefixFiles, ".txt"))
    doParallel::registerDoParallel(cluster)
  }

}

source("GenerateTestTreeAndData.R")

args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
  num_mpi_nodes <- as.integer(args[1])
} else {
  num_mpi_nodes <- 2
}

generatePCMModels <- function() {
  source("GeneratePCMModels.R")
}

modelTypes <- simulatedModels
argsMixedGaussian <- argsMixedGaussian_SimulatedModels
argsPCMParamLowerLimit <- list()
argsPCMParamUpperLimit <- list()

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)


tableFits <- NULL

#load("FitMappings_local_201812171_.RData")
#tableFits <- fitMappings$tableFits[, duplicated:=FALSE]

# fitsToClades <- PCMFitModelMappingsToClades(
#   values, tree.a, modelTypes = modelTypes,
#   generatePCMModelsFun = generatePCMModels,
#   metaIFun = PCMInfoCpp, positiveValueGuard = 1000,
#
#   tableFits = tableFits,
#
#   prefixFiles = prefixFiles,
#
#   minCladeSize = 30,
#
#   argsMixedGaussian = argsMixedGaussian,
#   argsPCMParamLowerLimit = argsPCMParamLowerLimit,
#   argsPCMParamUpperLimit = argsPCMParamUpperLimit,
#   argsConfigOptimAndMCMC = list(nCallsOptim = 2, genInitNumEvals = 1000, genInitVerbose = FALSE),
#
#   printFitVectorsToConsole = TRUE,
#   doParallel = TRUE,
#   verbose = TRUE)
#
# tableFits <- fitsToClades

fitMappings <- PCMFitModelMappings(
  values, tree.a, modelTypes = modelTypes,
  generatePCMModelsFun = generatePCMModels,
  metaIFun = PCMInfoCpp, positiveValueGuard = 1000,

  tableFits = tableFits,

  prefixFiles = prefixFiles,

  maxCladePartitionLevel = 10, maxNumNodesPerCladePartition = 1, minCladeSizes = 30,

  argsMixedGaussian = argsMixedGaussian,
  argsPCMParamLowerLimit = argsPCMParamLowerLimit,
  argsPCMParamUpperLimit = argsPCMParamUpperLimit,
  argsConfigOptimAndMCMC1 = list(nCallsOptim = 2, genInitNumEvals = 1000, genInitVerbose = FALSE),
  argsConfigOptimAndMCMC2 = list(nCallsOptim = 2, genInitNumEvals = 1000, genInitVerbose = FALSE),

  numJitterAllRegimeFits = 1000, numJitterRootRegimeFit = 1000,

  printFitVectorsToConsole = TRUE,
  doParallel = TRUE,
  verbose = TRUE)

save(fitMappings, file = paste0("FitMappings_", prefixFiles, ".RData"))

if(exists("cluster") && !is.null(cluster)) {
  parallel::stopCluster(cluster)
  # Don't forget to destroy the parallel cluster to avoid leaving zombie worker-processes.

  cluster <- NULL
}

# RETRIEVE THE BEST MODEL
#
# tableFits209 <- AddColumnFittedModel(fitMappings, fitMappings$tableFits[treeEDExpression=="E(tree,209)"][order(AIC)], TRUE)
#
# model1 <- tableFits209$fittedModel[[1]]
# model2 <- tableFits209$fittedModel[[2]]
# model3 <- tableFits209$fittedModel[[3]]
#
# bestModel <- AddColumnFittedModel(fitMappings, fitMappings$tableFits[treeEDExpression=="tree"][order(AIC)][1:2], TRUE)$fittedModel[[1]]
#
# plot <- PCMTreePlot(attr(bestModel, "tree")) + geom_nodelab(size = 2, color="black") +
#   theme(legend.position = "top")
# # View(tableFits[nobs==136][order(AIC)])
