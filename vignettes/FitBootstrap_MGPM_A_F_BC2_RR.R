# File: FitBootstrap_MGPM_A_F_BC2_RR.R
# Usage: R --vanilla --slave -f ../../FitBootstrap_MGPM_A_F_BC2_RR.R --args 1
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

# Creating the cluster for this PCMFit run.
# Uncomment the follwing code only if doMPI package is available and working.
# if(!exists("cluster") || is.null(cluster)) {
#   if(require(doMPI)) {
#     # using MPI cluster as distributed node cluster (possibly running on a
#     # cluster of multiple nodes)
#     # Get the number of cores. Assume this is run in a batch job.
#     p = strtoi(Sys.getenv('LSB_DJOB_NUMPROC'))
#     cluster <- startMPIcluster(count = p-1, verbose = TRUE)
#     doMPI::registerDoMPI(cluster)
#   } else {
#     # possibly running on personal computer without mpi installation
#     cluster <- parallel::makeCluster(
#       parallel::detectCores(logical = TRUE),
#       outfile = paste0("log_", prefixFiles, ".txt"))
#     doParallel::registerDoParallel(cluster)
#   }
#}

# This function is going to be executed on each worker node.
# Hence, it is the convenient place to define global settings,
# in particular:
# - set global option values,
# - call the function PCMGenerateModelTypes to generate the default PCM model types;
# - write custom S3 methods for PCMBase::PCMParentClasses and PCMBase::PCMSpecify
# for custom PCM model types (see corresponding man pages in the PCMBase package),
# - define S3 methods for the PCMBase::PCMParamUpperLimit and PCMParamLowerLimit in the
# global environment (see R-code below).
generatePCMModelsFunction <- function() {
  # make results reproducible but keep in mind that data_id is unique for this specific
  # worker node, and will be used to generate a unique parametric boostrap set of trait values.
  set.seed(data_id, kind = "Mersenne-Twister", normal.kind = "Inversion")

  PCMGenerateModelTypes()

  # An example DefineParameterLimits.R file can be found in
  # vignettes/DefineParameterLimits.R of the PCMFit package source-code.
  # Note that the path here is relative to the working directory of
  # the worker node R-process.
  source('../../DefineParameterLimits.R', local = FALSE)
}

# Inferred best model (
# replace this with the result from a previous call to  RetrieveBestFitScore(fitObject)$inferredModel)
bestModel <-
  RetrieveBestFitScore(PCMFitDemoObjects$fitMGPM_A_F_BC2_RR)$inferredModel
# Inferred tree
tree <- attr(bestModel, "tree")
# Simulate parameteric bootstrap trait values
X <- PCMSim(tree, bestModel, bestModel$X0)[, seq_len(PCMTreeNumTips(tree))]

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
