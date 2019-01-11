library(PCMFit)
library(PCMBase)

GeneratePCMModelTypes()

model2 <- do.call(
  MixedGaussian,
  c(list(
    k = 2,
    modelTypes = PCMFit::SurfaceOUTypeForMGPM(),
    mapping = c(a=1, b=1, c=1)),
    PCMFit::ArgsMPGM_SurfaceOU()))
