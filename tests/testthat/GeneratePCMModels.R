library(PCMBase)

listParameterizationsBM <- list(
  X0 = list(
    c("VectorParameter", "_Omitted")
  ),
  Sigma_x =
    list(
      c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal")
    ),

  Sigmae_x = list(
    # c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
    # c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
    # c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
    c("MatrixParameter", "_Omitted"))
)

listParameterizationsOU <- list(
  X0 = list(
    c("VectorParameter", "_Omitted")),

  H = list(
    c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
    c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
    c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
    c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),

    c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
    c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
    c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
    c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global")),

  Theta = list(
    c("VectorParameter")),

  Sigma_x = list(
    c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
    c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
    c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),

    c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
    c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
    c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global")),

  Sigmae_x = list(
    # c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
    # c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
    # c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
    c("MatrixParameter", "_Omitted"))
)

PCMGenerateParameterizations(structure(0.0, class="BM"), listParameterizations = listParameterizationsBM)
PCMGenerateParameterizations(structure(0.0, class="OU"), listParameterizations = listParameterizationsOU)


modelTypesMixedGaussian1 <- c(PCMModels("^BM__.*Omitted_Sigmae_x"), PCMModels("^OU__.*Schur.*Transformable_H.*Omitted_Sigmae_x"))
modelTypesMixedGaussian2 <- c(PCMModels("^BM__.*Omitted_Sigmae_x"), PCMModels("^OU__.*Schur.*Global_H.*Omitted_Sigmae_x"))
modelTypesMixedGaussian3 <- c(PCMModels("^OU__.*Schur.*ScalarDiagonal.*Global_H.*Omitted_Sigmae_x"))

simulatedModels <- c(
  # BM; independent traits
  "BM__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
  # BM; dependent traits
  "BM__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
  # OU; independent traits (diagonal H and diagonal Sigma)
  "OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
  # OU; dependtent traits (diagonal H and non-diagonal Sigma)
  "OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
  # OU; dependent traits (Symmetric H and non-diagonal Sigma)
  "OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
  # OU; dependent traits with causality (non-symmetric H and non-diagonal Sigma)
  "OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
)

argsMixedGaussian_SimulatedModels <- list(
  Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
                       description = "upper triangular Choleski factor of the non-phylogenetic variance-covariance")
)



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



inferredModel_SurfaceOU <- c(
  "OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x"
)
argsMixedGaussian_SurfaceOU <- list(
  H = structure(0.0,
                class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
                description = "adaptation rate matrix"),
  Sigma_x = structure(0.0,
                      class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
                      description = "unit-time variance parameter of the OU-process"),
  Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
                       description = "upper triangular Choleski factor of the non-phylogenetic variance-covariance")
)


inferredModel_ScalarOU <- c(
  "OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
)
argsMixedGaussian_ScalarOU <- list(
  H = structure(0.0,
                class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
                description = "adaptation rate matrix"),
  Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
                       description = "upper triangular Choleski factor of the non-phylogenetic variance-covariance")
)

inferredModel_MixedGaussian <- simulatedModels
argsMixedGaussian_MixedGaussian <- argsMixedGaussian_SimulatedModels
