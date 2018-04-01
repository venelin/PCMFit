library(ape)
library(testthat)
library(PCMBase)
library(PCMBaseCpp)
library(PCMFit)

set.seed(2)

# number of regimes
R <- 2
# number of traits
k <- 2

# number of tips
N <- 600

# rate mtrix of transition from one regime to another
Q <- matrix(c(-0.1, 0.1, 0.01, -0.01), R, R)
colnames(Q) <- rownames(Q) <- letters[1:R]


if(require(phytools)) {
  tree.a <- pbtree(n=N, scale=1, b = 1, d = 0.4)
  tree.ab <- phytools::sim.history(tree.a, Q, anc='a')
  tree.ab$edge.regime <- names(tree.ab$edge.length)

  # convert the simmap tree to a normal phylo object with singleton nodes at the
  # within-branch regime changes. The regimes are encoded as names of the edge.length
  # vector
  tree.ab.singles <- map.to.singleton(tree.ab)
  tree.ab.singles$edge.regime <- names(tree.ab.singles$edge.length)

  # insert a regime artificially at a random node
  PCMTreeSetRegimes(tree.ab.singles, nodes = 971, regimes = c("a", "b"))
  save(tree.ab.singles, file="tree.ab.singles.RData")
} else {
  load("tree.ab.singles.RData")
}
tree.ab.singles$edge.jump <- rep(0, nrow(tree.ab.singles$edge))

#PCMTreePlot(tree.ab.singles)


argsMRG <- list(X0 = list(default = rep(as.double(NA), 2), type = c("gvector", "fixed")),
                Sigmae_x = list(default = matrix(0, 2, 2), type = c("gmatrix", "fixed"),
                                description = "Fixed upper triangular Choleski factor of the variance-covariance matrix for the non-phylogenetic trait component"))
model <- do.call(MRG, c(list(k = 2, modelTypes = c("BM3", "OU4"), mapping = c(a = "BM3", b = "OU4")),
             argsMRG))

X0 <- c(5, 5)

a.Sigma_x <- rbind(
  c(1.6, 4),
  c(0, 2.4))

b.H <- rbind(
  c(2, 0),
  c(0, .6))
b.Theta <- c(4, 6)
b.Sigma_x <- rbind(
  c(1.6, 3),
  c(0, 0.3))


model$a$Sigma_x[,,1] <- a.Sigma_x

model$b$H[,,1] <- b.H
model$b$Theta[, 1] <- b.Theta
model$b$Sigma_x[,, 1] <- b.Sigma_x

traits.ab.123 <- PCMSim(tree.ab.singles, model, X0, verbose=TRUE)

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)
values <- traits.ab.123[, 1:length(tree.ab.singles$tip.label)]

likOrig <- PCMLik(values, tree.ab.singles, model)

options(PCMBase.Lmr.mode=21)

# cladeFits <- PCMFitModelMappings(
#   values, tree.ab.singles,
#   modelTypes = c("BM3", "OU4"),
#   fitToClades = TRUE,
#   printFitVectorsToConsole = TRUE,
#   doParallel = FALSE,
#   nCallsOptim = 100,
#   genInitVerbose = FALSE, verbose = TRUE, metaI = PCMInfoCpp,
#   X0 = list(default = rep(as.double(NA), 2), type = c("gvector", "fixed")),
#   Sigmae_x = list(default = matrix(0, 2, 2), type = c("gmatrix", "fixed")))

load("cladeFits.RData")

# if(!exists("cluster") || is.null(cluster)) {
#   cluster <- parallel::makeCluster(parallel::detectCores(logical = TRUE), outfile = "log.txt")
#   doParallel::registerDoParallel(cluster)
# }
load("select.RData")

# select <- PCMFitModelMappings(
#   values, tree.ab.singles,
#   modelTypes = c("BM3", "OU4"),
#   printFitVectorsToConsole = TRUE,
#   doParallel = FALSE,
#   nCallsOptim = 100,
#   genInitVerbose = FALSE, verbose = TRUE, metaI = PCMInfoCpp,
#   X0 = list(default = rep(as.double(NA), 2), type = c("gvector", "fixed")),
#   Sigmae_x = list(default = matrix(0, 2, 2), type = c("gmatrix", "fixed")))


select2 <- garbage<-PCMFitModelMappings(
  X = values, tree = tree.ab.singles, modelTypes = c("BM3", "OU4"), metaIFun = PCMInfoCpp,


  cladeFits = cladeFits,

  numJitterAllCladeFits = 1000, numJitterRootCladeFit = 1000,
  argsMRG = argsMRG,
  argsConfigOptimAndMCMC = list(nCallsOptim = 10, genInitNumEvals = 10, genInitVerbose = FALSE),

  printFitVectorsToConsole = TRUE,
  doParallel = FALSE,
  verbose = TRUE
  )

# select3 <- PCMFitModelMappings(
#   values, tree.ab.singles,
#   modelTypes = c("BM3", "OU4"),
#   printFitVectorsToConsole = TRUE,
#   cladeFits = cladeFits,
#   doParallel = FALSE,
#   nCallsOptim = 100,
#   genInitVerbose = FALSE, verbose = TRUE, metaI = PCMInfoCpp,
#   X0 = list(default = rep(as.double(NA), 2), type = c("gvector", "fixed")),
#   Sigmae_x = list(default = matrix(0, 2, 2), type = c("gmatrix", "fixed")))


modelBestAIC <- PCMLoadMRGFromFitVector(fitVector = garbage[[2]], k = 2, modelTypes = c("BM3", "OU4"),
                                        X0 = list(default = rep(as.double(NA), 2), type = c("gvector", "fixed")),
                                        Sigmae_x = list(default = matrix(0, 2, 2), type = c("gmatrix", "fixed")))
#
# modelClade971 <- PCMLoadMRGFromFitVector(cladeFits[[2]][[2]],
#                                          k = 2,
#                                          modelTypes = c("BM3", "OU4"),
#                                          X0 = list(default = rep(as.double(NA), 2), type = c("gvector", "fixed")),
#                                          Sigmae_x = list(default = matrix(0, 2, 2), type = c("gmatrix", "fixed")))
#
# treeSplit971 <- PCMTreeSplitAtNode(tree.ab.singles, 971, X=values)
# clade971 <- treeSplit971$clade
# PCMTreeSetDefaultRegime(clade971, 1)
# X971 <- treeSplit971$Xclade

# fit971_OU4 <- PCMFit(X971, clade971,
#                      model = PCMLoadMRGFromFitVector(cladeFits[[2]][[2]], k=2, modelTypes = c("BM3", "OU4"),
#                                                      X0 = list(default = rep(as.double(NA), 2), type = c("gvector", "fixed")),
#                                                      Sigmae_x = list(default = matrix(0, 2, 2), type = c("gmatrix", "fixed"))),
#                      nCallsOptim = 200, genInitVerbose = TRUE, verbose = TRUE,
#                      metaI = PCMInfoCpp)

# modelFit971_OU4 <- fit971_OU4$modelOptim
# model971_OU4_orig <- modelFit971_OU4
# PCMLik(X971, clade971, modelFit971_OU4)
#
# model971_OU4_orig[["1"]]$H[,,1]<-b.H
# model971_OU4_orig[["1"]]$Theta[,1]<-b.Theta
# model971_OU4_orig[["1"]]$Sigma_x[,,1]<-b.Sigma_x
# PCMLik(X971, clade971, model971_OU4_orig)


# Don't forget to destroy the parallel cluster to avoid leaving zombie worker-processes.
# parallel::stopCluster(cluster)
# cluster <- NULL

for(res in select) {
  cat(toString(unname(res)), "\n")
}

# test_that("ML value is better than original",
#           expect_true(logLik(fit2) < 0 && logLik(fit2) > likOrig))
