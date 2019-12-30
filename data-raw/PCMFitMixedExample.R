# running the example on a cluster:
# bsub -M 100000 -n 24 -W 23:59 -R ib sh R --vanilla --slave -f ../PCMFitMixedExample.R
# running locally using a shell:
# sh R --vanilla --slave -f ../PCMFitMixedExample.R

library(PCMBase)
library(PCMBaseCpp)
library(PCMFit)
library(data.table)
# other needed packages, e.g. ape, data.table etc...

# A character string used in filenames for a model inference on a given data:
prefixFiles = paste0("MGPM_A_F_BC2_RR")

# creating the cluster for this PCMFit run:
if(!exists("cluster") || is.null(cluster)) {
  if(require(doMPI)) {
    # using MPI cluster as distributed node cluster (possibly running on a
    # cluster of multiple nodes)
    # Get the number of cores. Assume this is run in a batch job.
    p = strtoi(Sys.getenv('LSB_DJOB_NUMPROC'))
    cluster <- startMPIcluster(count = p-1, verbose = TRUE)
    doMPI::registerDoMPI(cluster)
  } else if(require(parallle)) {
    # possibly running on personal computer without mpi installation
    cluster <- parallel::makeCluster(
      parallel::detectCores(logical = TRUE),
      outfile = paste0("log_", prefixFiles, ".txt"))
    doParallel::registerDoParallel(cluster)
  } else {
    warning("Could not create a cluster. Running serially.")
  }
}

# This function is going to be executed on each worker node.
generatePCMModelsFunction <- function() {
  # make results reproducible
  set.seed(4, kind = "Mersenne-Twister", normal.kind = "Inversion")

  PCMGenerateModelTypes()
  fileName <- '../DefineParameterLimits.R'
  codeDefineLimits <- readChar(fileName, file.info(fileName)$size)
  eval(parse(text = codeDefineLimits), .GlobalEnv)
}

tree <- PCMTree(PCMFitDemoObjects$dtSimulated$tree[[1]])
X <- PCMFitDemoObjects$dtSimulated$X[[1]][, seq_len(PCMTreeNumTips(tree))]

currentResultFile <- paste0("Current_", prefixFiles, ".RData")
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

bestFit <- RetrieveBestFitScore(fitMGPM_A_F_BC2_RR)

modelTrue <- PCMFitDemoObjects$dtSimulated$model[[1]]

# We specify the tree and trait values for the true model in order to easily
# calculate parameter count likelihood and AIC for it:
attr(modelTrue, "tree") <- PCMFitDemoObjects$dtSimulated$treeWithRegimes[[1]]
attr(modelTrue, "X") <- X
attr(modelTrue, "SE") <- X * 0.0

listModels <- list(
  RetrieveBestModel(PCMFitDemoObjects$fitBM),
  RetrieveBestModel(PCMFitDemoObjects$fitOU),
  RetrieveBestModel(PCMFitDemoObjects$fitMGPMTrueTypeMapping),
  RetrieveBestModel(PCMFitDemoObjects$fitMGPMTrueTypeMappingCheat),
  modelTrue,
  bestFit$inferredModel)

dtSummary <- data.table(
  model = c(
    "Global BM",
    "Global OU",
    "True MGPM, unknown parameters",
    "True MGPM, known true parameters",
    "True MGPM, true parameters",
    "Inferred MGPM, unknown shift-points and parameters"),
  p = sapply(listModels, PCMParamCount),
  logLik = sapply(listModels, logLik),
  AIC = sapply(listModels, AIC))

save(fitMGPM_A_F_BC2_RR, bestFit, listModels, dtSummary,
     file = "ResultPCMFitExample.RData")
