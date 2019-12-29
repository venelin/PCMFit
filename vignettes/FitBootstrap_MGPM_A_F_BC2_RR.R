library(PCMBase)
library(PCMBaseCpp)
library(PCMFit)
library(data.table)

# extract dataset identifier and possibly other parameters from the command line:
args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
  data_id <- as.integer(args[1])
} else {
  data_id <- 1L
}

# A character string used in filenames for a model inference on a given data:
prefixFiles = paste0("MGPM_A_F_BC2_RR_BSID_", data_id)

# creating the cluster for this PCMFit run:
if(!exists("cluster") || is.null(cluster)) {
  if(require(doMPI)) {
    # using MPI cluster as distributed node cluster (possibly running on a
    # cluster of multiple nodes)
    # Get the number of cores. Assume this is run in a batch job.
    p = strtoi(Sys.getenv('LSB_DJOB_NUMPROC'))
    cluster <- startMPIcluster(count = p-1, verbose = TRUE)
    doMPI::registerDoMPI(cluster)
  } else {
    # possibly running on personal computer without mpi installation
    cluster <- parallel::makeCluster(
      parallel::detectCores(logical = TRUE),
      outfile = paste0("log_", prefixFiles, ".txt"))
    doParallel::registerDoParallel(cluster)
  }
}

# This function is going to be executed on each worker node.
generatePCMModelsFunction <- function() {
  # make results reproducible
  set.seed(4, kind = "Mersenne-Twister", normal.kind = "Inversion")

  PCMGenerateModelTypes()
  fileName <- '../../DefineParameterLimits.R'
  codeDefineLimits <- readChar(fileName, file.info(fileName)$size)
  eval(parse(text = codeDefineLimits), .GlobalEnv)
}

tree <- PCMTree(PCMFitDemoObjects$dtSimulated$tree[[1]])
X <- PCMFitDemoObjects$valuesBootstrapMGPM_A_F_BC2_RR$X[[data_id]][
  , seq_len(PCMTreeNumTips(tree))]

currentResultFile <- paste0("CurrentResults_fits_", prefixFiles, ".RData")
if(file.exists(currentResultFile)) {
  load(currentResultFile)
  tableFitsPrev <- listResults$tableFits
} else {
  tableFitsPrev <- NULL
}

fitMGPM_A_F_BC2_RR <- PCMFitMixed(
  X = X, tree = tree, metaIFun = PCMInfoCpp,
  generatePCMModelsFun = generatePCMModelsFunction,
  maxNumRoundRobins = 2, maxNumPartitionsInRoundRobins = 2,
  tableFitsPrev = tableFitsPrev,
  prefixFiles = prefixFiles,
  doParallel = TRUE)

save(fitMGPM_A_F_BC2_RR, file = paste0("Result_", prefixFiles, ".RData"))
