modelTypes <- PCMModels("^(BM|OU).*__noX0.*__noSigmae_x")

initMapping = c("BM__noX0__noSigmae_x", "BM__noX0__noSigmae_x", "BM__noX0__noSigmae_x__posdiagSigma_x",
                "OU__noX0__noSigmae_x__posdiagH__posdiagSigma_x")

allowedModelTypesIndices <- list(NULL, c(1,3), c(2,4), 5)

it <- PCMIteratorMapping2(mapping = initMapping, modelTypes = modelTypes, allowedModelTypesIndices = allowedModelTypesIndices)

library(iterators)
for(i in 1:20) cat(toString(nextElem(it)), "\n")
