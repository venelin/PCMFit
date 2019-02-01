library(PCMBase)




inferredModel_BM <- c(
  "BM__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
)

argsMixedGaussian_BM <- list(
  Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
                       description = "upper triangular Choleski factor of the non-phylogenetic variance-covariance")
)




inferredModel_OU <- c(
  "OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
)

argsMixedGaussian_OU <- list(
  Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
                       description = "upper triangular Choleski factor of the non-phylogenetic variance-covariance")
)






inferredModel_MixedGaussian <- simulatedModels
argsMixedGaussian_MixedGaussian <- argsMixedGaussian_SimulatedModels
