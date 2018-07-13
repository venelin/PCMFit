library(PCMBase)
library(PCMFit)

listParameterizationsBM <- list(
  X0 = list(c("VectorParameter", "_Omitted")),
  Sigma_x = list(c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
                 c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
                 c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal")),

  Sigmae_x = list(c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
                  c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
                  c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
                  c("MatrixParameter", "_Omitted"))
)

listParameterizationsOU <- list(
  X0 = list(c("VectorParameter", "_Omitted")),
  H = list(c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
           c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
           c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
           c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),

           c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
           c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
           c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
           c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global")),

  Theta = list(c("VectorParameter")),

  Sigma_x = list(c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
                 c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
                 c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal")),

  Sigmae_x = list(c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
                  c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
                  c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
                  c("MatrixParameter", "_Omitted"))
)

PCMGenerateParameterizations(structure(0.0, class="BM"), listParameterizations = listParameterizationsBM)
PCMGenerateParameterizations(structure(0.0, class="OU"), listParameterizations = listParameterizationsOU)

modelTypes <- PCMModels()


initMapping = c("BM__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigmae_x",
                "BM__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
                "BM__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__ScalarDiagonal_WithNonNegativeDiagonal_Sigmae_x",
                "OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__ScalarDiagonal_WithNonNegativeDiagonal_Sigmae_x")

allowedModelTypesIndices <- list(NULL, c(1:30), c(2:50), c(2:60))

it <- PCMIteratorMapping2(mapping = initMapping, modelTypes = modelTypes, allowedModelTypesIndices = allowedModelTypesIndices)
library(iterators)
res <- try(
for(i in 1:20) cat(toString(nextElem(it)), "\n")
)
