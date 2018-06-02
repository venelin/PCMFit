library(data.table)

if(!exists("cluster") || is.null(cluster)) {
  if(require(doMPI)) {
    # using MPI cluster as distributed node cluster (possibly running on a cluster)
    # Get the number of cores. Assume this is run in a batch job.
    p = strtoi(Sys.getenv('LSB_DJOB_NUMPROC'))
    cluster <- startMPIcluster(count = p-1, verbose = TRUE)
    doMPI::registerDoMPI(cluster)
  } else {
    cluster <- parallel::makeCluster(parallel::detectCores(logical = TRUE), outfile = "log27.txt")
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

modelTypes <-  c(
  "BM__NoX0__NoSigmae_x",
  #"BM__NoX0__NoSigmae_x__DiagPosdiagSigma_x",
  "OU__NoX0__NoSigmae_x__PosdiagH"#,
  #"OU__NoX0__NoSigmae_x__SymmetricPosdiagH",
  #"OU__NoX0__NoSigmae_x__DiagPosdiagH",
  #"OU__NoX0__NoSigmae_x__DiagPosdiagH__DiagPosdiagSigma_x"#,
  #"OU__NoX0__NoSigmae_x__PosdiagH__DiagPosdiagSigma_x"#,
  # "OU__NoX0__NoSigmae_x__SymmetricPosdiagH__DiagPosdiagSigma_x",
  # "OU__NoX0__NoSigmae_x__DiagPosdiagH__DiagPosdiagSigma_x",
  # "DOU__NoX0__NoSigmae_x__ZeroH1__PosdiagH2",
  # "DOU__NoX0__NoSigmae_x__ZeroH1__SymmetricPosdiagH2",
  # "DOU__NoX0__NoSigmae_x__ZeroH1__DiagPosdiagH2",
  # "DOU__NoX0__NoSigmae_x__ZeroH1__PosdiagH2__DiagPosdiagSigma_x",
  # "DOU__NoX0__NoSigmae_x__ZeroH1__SymmetricPosdiagH2__DiagPosdiagSigma_x",
  # "DOU__NoX0__NoSigmae_x__ZeroH1__DiagPosdiagH2__DiagPosdiagSigma_x",
  # "DOU__NoX0__NoSigmae_x__DiagPosdiagH1__PosdiagH2",
  # "DOU__NoX0__NoSigmae_x__DiagPosdiagH1__SymmetricPosdiagH2",
  # "DOU__NoX0__NoSigmae_x__DiagPosdiagH1__DiagPosdiagH2",
  # "DOU__NoX0__NoSigmae_x__DiagPosdiagH1__PosdiagH2__DiagPosdiagSigma_x",
  # "DOU__NoX0__NoSigmae_x__DiagPosdiagH1__SymmetricPosdiagH2__DiagPosdiagSigma_x",
  # "DOU__NoX0__NoSigmae_x__DiagPosdiagH1__DiagPosdiagH2__DiagPosdiagSigma_x"
  )

argsMRG <- list(X0 = list(default = rep(as.double(0), 2), type = c("gvector", "full")),
                Sigmae_x = list(default = matrix(0, 2, 2), type = c("gmatrix", "fixed"),
                                description = "Fixed upper triangular Choleski factor of the variance-covariance matrix for the non-phylogenetic trait component"))

argsPCMLowerBound <- list(lowerBoundValue = -10, lowerBoundValuePositiveDiag = 0)
argsPCMUpperBound <- list(upperBoundValue = 10, upperBoundValuePositiveDiag = 10)

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)


tableFits <- NULL
load("fitMappings_201805194.RData")
tableFits <- fitMappings$tableFits
# tableFits <- rbindlist(lapply(list.files(".", pattern="fits_201805192_.*.RData"), function(f)  {
#   load(f)
#   data
# }))

print(tableFits)
tableFits[, duplicated:=FALSE]

if(nrow(tableFits) == 0) {
  tableFits <- NULL
}


fitMappings <- PCMFitModelMappings(
  values, tree.ab.singles, modelTypes = modelTypes,
  metaIFun = PCMInfoCpp, positiveValueGuard = 1000,

  tableFits = tableFits,

  prefixFiles = "fits_201805195_",
  maxCladePartitionLevel = 3, minCladeSizes = c(38, 20, 15),

  argsMRG = argsMRG,
  argsPCMLowerBound = argsPCMLowerBound,
  argsPCMUpperBound = argsPCMUpperBound,
  argsConfigOptimAndMCMC1 = list(nCallsOptim = 25, genInitNumEvals = 1000, genInitVerbose = FALSE),
  argsConfigOptimAndMCMC2 = list(nCallsOptim = 4, genInitNumEvals = 1000, genInitVerbose = FALSE),

  numJitterAllRegimeFits = 1000, numJitterRootRegimeFit = 1000,

  printFitVectorsToConsole = TRUE,
  doParallel = TRUE,
  verbose = TRUE)

save(fitMappings, file = "fitMappings_201805195.RData")

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
