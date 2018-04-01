library(ape)
library(testthat)
library(PCMBase)
library(PCMFit)
library(abind)

set.seed(2)

# number of regimes
R <- 2
# number of traits
k <- 2

# number of tips
N <- 80

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
} else {
  tree.a <- rtree(N)
  tree.ab <- tree.a
  tree.ab$edge.regime <- sample(c("a", "b"), size = length(tree.ab$edge.length), replace = TRUE)
  tree.ab.singles <- tree.ab
}
tree.ab.singles$edge.jump <- rep(0, nrow(tree.ab.singles$edge))

#PCMPlotTree(tree.ab.singles)

model <- MRG(k = 2, modelTypes = c("BM3", "OU3"), mapping = c(a = "BM3", b = "OU3"),
             X0 = list(default = rep(as.double(NA), 2), type = c("gvector", "fixed")),
             Sigmae_x = list(default = matrix(0, 2, 2), type = c("gmatrix", "fixed"),
                           description = "Fixed upper triangular Choleski factor of the variance-covariance matrix for the non-phylogenetic trait component"))

X0 <- c(5, 5)

a.Sigma_x <- rbind(
  c(1.6, 4),
  c(0, 2.4))

b.H <- rbind(
  c(2, 2),
  c(0, .6))
b.Theta <- c(10, 6)
b.Sigma_x <- rbind(
  c(1.6, 3),
  c(0, 0.3))


model$a$Sigma_x[,,1] <- a.Sigma_x

model$b$H[,,1] <- b.H
model$b$Theta[, 1] <- b.Theta
model$b$Sigma_x[,, 1] <- b.Sigma_x

traits.ab.123 <- PCMSim(tree.ab.singles, model, X0, verbose=TRUE)

options(PCMBase.Value.NA = -1e20)
values <- traits.ab.123[, 1:length(tree.ab.singles$tip.label)]

likOrig <- PCMLik(values, tree.ab.singles, model)

if(require(PCMBaseCpp)) {
  metaInfo <- PCMInfoCpp(values, tree.ab.singles, model)
} else {
  metaInfo <- PCMInfo(values, tree.ab.singles, model)
}

likFun2 <- PCMCreateLikelihood(values, tree.ab.singles, model, metaInfo)

fit2 <- PCMFit(values, tree.ab.singles, model, metaI = metaInfo,
               argsConfigOptimAndMCMC = list(nCallsOptim = 20, genInitVerbose = TRUE),
               verbose = TRUE)

test_that("ML value is better than original",
          expect_true(logLik(fit2) < 0 && logLik(fit2) > likOrig))
