library(ape)
library(testthat)
library(PCMBase)
library(PCMStep)
library(abind)

set.seed(2)

# number of regimes
R <- 2
# number of traits
k <- 2

# number of tips
N <- 100

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

PCMPlotTree(tree.ab.singles)

model <- MRG(k = 2, models = c("BM3", "OU3"), mapping = c(a = "BM3", b = "OU3"),
             X0 = list(default = rep(as.double(NA), 2), type = c("gvector", "fixed")),
             Sigmae = list(default = matrix(0, 2, 2), type = c("gmatrix", "fixed"),
                           description = "variance-covariance matrix for the non-phylogenetic trait component"))

X0 <- c(5, 5)

a.Sigma <- rbind(
  c(1.6, 0),
  c(0, 2.4))

b.H <- rbind(
  c(2, .1),
  c(.1, .6))
b.Theta <- c(10, 6)
b.Sigma <- rbind(
  c(1.6, .3),
  c(.3, 0.3))

Sigmae <- rbind(
  c(.2, 0),
  c(0, .3))


model$a$Sigma[,,1] <- a.Sigma

model$b$H[,,1] <- b.H
model$b$Theta[, 1] <- b.Theta
model$b$Sigma[,, 1] <- b.Sigma

#model$Sigmae <- Sigmae

traits.ab.123 <- PCMSim(tree.ab.singles, model, X0, verbose=TRUE)

options(PCMBase.Value.NA = -1e20)
values <- traits.ab.123[, 1:length(tree.ab.singles$tip.label)]


PCMLik(values, tree.ab.singles, model)

likFun <- PCMCreateLikelihood(values, tree.ab.singles, model)

if(require(PCMBaseCpp)) {
  metaInfo <- PCMInfoCpp(values, tree.ab.singles, model)
} else {
  metaInfo <- PCMInfo(values, tree.ab.singles, model)
}

likFun2 <- PCMCreateLikelihood(values, tree.ab.singles, model, metaInfo)

#fit <- PCMFit(values, tree.ab.singles, model)
library(OptimMCMC)
# listParInitOptim = genInitStates(rep(-10, 12), rep(10, 12), 10, likFun2, minValue = -1e20, verbose = TRUE)

fit2 <- PCMFit(values, tree.ab.singles, model, metaI = metaInfo, nCallsOptim = 500,
               minGenInitValue = -1e20, maxGenInitIter = 5000, verboseGenInit = TRUE, verbose = TRUE)

model2 <- model
PCMSetOrGetVecParams(model2, fit2$Optim$par)
print(model2)

# toString(likFun2(fit2$Optim$par))
#

if(FALSE) {


  par <- c(9.20387191445074, -3.36841497851032, 9.07657202847133, 5.73227273365591, -8.11659779441043, 6.81959596230008, -1.7395069588348, -0.39438765743306, -1.78697063300318, 8.58882728994634, -1.50839757416005, 8.06732489221055)

  likFun(par)
  likFun2(par)

  model2 <- model

  PCMSetOrGetVecParams(model2, par)
  model2

  par2 <- par
  PCMSetOrGetVecParams(model, par2, set=FALSE)

  Lmr<-PCMLmr(X = values, tree = tree.ab.singles, model = model2)
  Lmr2<-PCMLmr(X = values, tree = tree.ab.singles, model = model2, metaInfo)

  metaInfoR <- PCMInfo(values, tree.ab.singles, model2)
  condOU <- PCMCond(tree = tree.ab.singles, model = model2, r = 2)
  edgeIndex <- which(tree.ab.singles$edge[, 2]==8)
  time <- tree.ab.singles$edge.length[edgeIndex]

  V <- condOU$V(time, edgeIndex, metaInfoR)


  lstR <- PCMBase:::PCMPLambdaP_1(model2$b$H[,,1])
  Lambda_ij <- PCMBase:::PCMPairSums(lstR$lambda)
  fLambda_ij <- PCMBase:::PCMPExpxMeanExp(Lambda_ij)(time)

  P_1SigmaP_1_t <- lstR$P_1 %*% model2$b$Sigma[,,1] %*% t(lstR$P_1)

  lstR$P %*% (fLambda_ij * P_1SigmaP_1_t) %*% t(lstR$P)

  lst <- PCMBaseCpp:::VOU(model2$b$H[,,1], model2$b$Sigma[,,1], time, 1e-6, 1e-6)
  lst$P <- lst$P_real + lst$P_imag*1i
  lst$P_1 <- lst$P_1_real + lst$P_1_imag*1i
  lst$fLambda_ij  <- lst$fLambda_ij_real + lst$fLambda_ij_imag*1i
  lst$P_1SigmaP_1_t <- lst$P_1SigmaP_1_t_real + lst$P_1SigmaP_1_t_imag*1i


  lst$Sigma - model2$b$Sigma[,,1]
  lst$P - lstR$P
  lst$P_1 - lstR$P_1
  lst$fLambda_ij - fLambda_ij
  lst$P_1SigmaP_1_t - P_1SigmaP_1_t

  lst$P_1_real %*% model2$b$Sigma[,,1] %*% t(lst$P_1_real)

  lst$P_1 %*% (model2$b$Sigma[,,1] %*% t(lst$P_1))

  lst$P %*% (lst$fLambda_ij * lst$P_1SigmaP_1_t) %*% t(lst$P)
  lst$P %*% (lst$fLambda_ij * P_1SigmaP_1_t) %*% t(lst$P)
}
