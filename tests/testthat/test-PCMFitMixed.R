library(data.table)
library(PCMFit)
prefixFiles <- "local_20190111_"

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
  PCMFit::GeneratePCMModelTypes()
}

modelTypes <- PCMFit::MGPMDefaultModelTypes()
argsMixedGaussian <- PCMFit::Args_MixedGaussian_MGPMDefaultModelTypes()
argsPCMParamLowerLimit <- list()
argsPCMParamUpperLimit <- list()

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)


tableFits <- NULL

load("FitMappings_local_201901106_.RData")
tableFits <- fitMappings$tableFits[, duplicated:=FALSE]
#setnames(tableFits, "aic", "score")

fitMappings <- PCMFitMixed(
  values, tree.a, modelTypes = modelTypes,
  generatePCMModelsFun = generatePCMModels,
  metaIFun = PCMInfoCpp, positiveValueGuard = 1000,

  tableFits = tableFits,

  #listPartitions = list(c(101, 108, 105, 115)),
  listPartitions = "all",

  minCladeSizes = 20,

  maxCladePartitionLevel = 1, maxNumNodesPerCladePartition = Inf,
  #maxCladePartitionLevel = 100, maxNumNodesPerCladePartition = 1,


  #listAllowedModelTypesIndices = "all",
  listAllowedModelTypesIndices = "best-clade",

  argsMixedGaussian = argsMixedGaussian,

  argsConfigOptim1 = defaultArgsConfigOptim(
    numCallsOptim = 2,
    argsPCMParamLowerLimit = argsPCMParamLowerLimit,
    argsPCMParamUpperLimit = argsPCMParamUpperLimit),
  argsConfigOptim2 = defaultArgsConfigOptim(
    numCallsOptim = 2,
    argsPCMParamLowerLimit = argsPCMParamLowerLimit,
    argsPCMParamUpperLimit = argsPCMParamUpperLimit),

  doParallel = TRUE,

  prefixFiles = prefixFiles,
  saveTempWorkerResults = TRUE,
  printFitVectorsToConsole = TRUE,
  verbose = TRUE,
  debug = FALSE)

save(fitMappings, file = paste0("FitMappings_", prefixFiles, ".RData"))
if(exists("cluster") && !is.null(cluster)) {
  parallel::stopCluster(cluster)
  # Don't forget to destroy the parallel cluster to avoid leaving zombie worker-processes.

  cluster <- NULL
}

