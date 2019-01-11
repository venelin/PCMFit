library(ape)
library(testthat)
library(PCMBase)

library(abind)
library(data.table)

PCMFit::GeneratePCMModelTypes()

set.seed(2)


# number of regimes
R <- 2
# number of traits
k <- 2

# number of tips
N <- 100

# rate mtrix of transition from one regime to another
# Q <- matrix(c(-0.8, 0.8, 0.01, -0.01), R, R)
# colnames(Q) <- rownames(Q) <- letters[1:R]


tree.a <- rcoal(n=N)

PCMTreeSetDefaultRegime(tree.a, "a")
PCMTreeSetRegimes(tree.a, c(104, 102), regimes = c("a", "b", "c"))
PCMTreeSetLabels(tree.a)


tree.a$edge.jump <- rep(0, nrow(tree.a$edge))

#PCMTreePlot(tree.a) + geom_nodelab()


model1 <- do.call(
  MixedGaussian,
  c(list(
    k = 2,
    modelTypes = PCMFit::MGPMDefaultModelTypes(),
    mapping = c(a=1, b=4, c=5)),
    PCMFit::Args_MixedGaussian_MGPMDefaultModelTypes()))



vecModelRandom <- round(PCMParamRandomVecParams(model1), 1)
modelRandom <- model1
PCMParamLoadOrStore(modelRandom, vecModelRandom, offset = 0, load=TRUE)
model1 <- modelRandom

#model1$c$H[1,2,1] <- 0

traits <- PCMSim(tree.a, model1, X0 = model1$X0, verbose=TRUE)

options(PCMBase.Value.NA = -1e20)
values <- traits[, 1:length(tree.a$tip.label)]

likR <- PCMLik(values, tree.a, model1)

if(require(PCMBaseCpp)) {
  metaInfo <- PCMInfoCpp(values, tree.a, model1)
} else {
  metaInfo <- PCMInfo(values, tree.a, model1)
}

likC <- PCMLik(values, tree.a, model1, metaI = PCMInfoCpp(values, tree.a, model1))

likRT <- PCMLik(values, tree.a, PCMApplyTransformation(model1))
likCT <- PCMLik(values, tree.a, PCMApplyTransformation(model1), metaI = PCMInfoCpp(values, tree.a, PCMApplyTransformation(model1)))

#likFun2 <- PCMCreateLikelihood(values, tree.a, model1, metaInfo)

#print(likFun2(p = PCMParamGetShortVector(model1)))
print(likC)
print(likR)

#
# tree.a$edge[tree.a$edge[,1]==390,]
#
# LmrR<-PCMLmr(values, tree.a, PCMApplyTransformation(model1))
# LmrC<-PCMLmr(values, tree.a, PCMApplyTransformation(model1), metaI = PCMInfoCpp(values, tree.a, PCMApplyTransformation(model1)))
