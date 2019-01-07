#' @export
defaultArgsConfigOptim <- function(
  matParInit = NULL,
  argsPCMParamLowerLimit = NULL,
  argsPCMParamUpperLimit = NULL,
  numRunifInitVecParams = 1000,
  numGuessInitVecParams = 100, varyThetaInGuessInitVecParams = FALSE,
  numJitterRootRegimeFit = 100, sdJitterRootRegimeFit = 0.5,
  numJitterAllRegimeFits = 100, sdJitterAllRegimeFits = 0.5,
  numCallsOptim = 10,
  control = list(fnscale = -1) ) {

  as.list(environment())
}

#' @export
ListDefaultParameterizations <- function(o) {
  UseMethod("ListDefaultParameterizations", o)
}

#' @export
ListDefaultParameterizations.BM <- function(o) {
  list(
    X0 = list(
      c("VectorParameter", "_Global"),
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
}

#' @export
ListDefaultParameterizations.OU <- function(o) {
  list(
    X0 = list(
      c("VectorParameter", "_Global"),
      c("VectorParameter", "_Omitted")
    ),

    H = list(
      c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable"),
      c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),
      c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable"),
      c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable"),

      c("MatrixParameter", "_Schur", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
      c("MatrixParameter", "_Schur", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
      c("MatrixParameter", "_Schur", "_Diagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
      c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global")
    ),

    Theta = list(
      c("VectorParameter")),

    Sigma_x = list(
      c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),

      c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
      c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
      c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global")
    ),

    Sigmae_x = list(
      # c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
      # c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
      # c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_Omitted"))
  )
}

#' @export
GeneratePCMModelTypes <- function(baseTypes = c("BM", "OU")) {
  for(bt in baseTypes) {
    o <- structure(0.0, class=bt)
    PCMGenerateParameterizations(
      o,
      listParameterizations = ListDefaultParameterizations(o))
  }
}


#' @export
DefaultModelTypes <- function() {
  c(
    # BM; independent traits
    A = "BM__Global_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    # BM; dependent traits
    B = "BM__Global_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    # OU; independent traits (diagonal H and diagonal Sigma)
    C = "OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    # OU; dependtent traits (diagonal H and non-diagonal Sigma)
    D = "OU__Global_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    # OU; dependent traits (Symmetric H and non-diagonal Sigma)
    E = "OU__Global_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    # OU; dependent traits with causality (non-symmetric H and non-diagonal Sigma)
    F = "OU__Global_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
  )
}

#' @export
DefaultModelTypesForMGPM <- function() {
  c(
    # BM; independent traits
    A = "BM__Omitted_X0__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    # BM; dependent traits
    B = "BM__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    # OU; independent traits (diagonal H and diagonal Sigma)
    C = "OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__Diagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    # OU; dependtent traits (diagonal H and non-diagonal Sigma)
    D = "OU__Omitted_X0__Schur_Diagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    # OU; dependent traits (Symmetric H and non-diagonal Sigma)
    E = "OU__Omitted_X0__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    # OU; dependent traits with causality (non-symmetric H and non-diagonal Sigma)
    F = "OU__Omitted_X0__Schur_WithNonNegativeDiagonal_Transformable_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
  )
}
#' @export
ArgsMGPM_DefaultModelTypes <- function() {
  list(
    Sigmae_x = structure(
      0.0, class = c("MatrixParameter", "_Omitted"),
      description = "upper triangular Choleski factor of the non-phylogenetic variance-covariance")
  )
}


#' @export
ScalarOUTypeForMGPM <- function() {
  c(
    "OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"
  )
}

#' @export
ArgsMGPM_ScalarOU <- function() {
  list(
    H = structure(
      0.0,
      class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
      description = "adaptation rate matrix"),
    Sigmae_x = structure(0.0, class = c("MatrixParameter", "_Omitted"),
                         description = "upper triangular Choleski factor of the non-phylogenetic variance-covariance")
  )
}

#' @export
SurfaceOUTypeForMGPM <- c(
  "OU__Omitted_X0__Schur_ScalarDiagonal_WithNonNegativeDiagonal_Transformable_Global_H__Theta__Diagonal_WithNonNegativeDiagonal_Global_Sigma_x__Omitted_Sigmae_x"
)
#' @export
ArgsMPGM_SurfaceOU <- list(
  H = structure(
    0.0,
    class = c("MatrixParameter", "_Schur", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Transformable", "_Global"),
    description = "adaptation rate matrix"),
  Sigma_x = structure(
    0.0,
    class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
    description = "unit-time variance parameter of the OU-process"),
  Sigmae_x = structure(
    0.0, class = c("MatrixParameter", "_Omitted"),
    description = "upper triangular Choleski factor of the non-phylogenetic variance-covariance")
)

#
# modelTypesMixedGaussian1 <- c(PCMModels("^BM__.*Omitted_Sigmae_x"), PCMModels("^OU__.*Schur.*Transformable_H.*Omitted_Sigmae_x"))
# modelTypesMixedGaussian2 <- c(PCMModels("^BM__.*Omitted_Sigmae_x"), PCMModels("^OU__.*Schur.*Global_H.*Omitted_Sigmae_x"))
# modelTypesMixedGaussian3 <- c(PCMModels("^OU__.*Schur.*ScalarDiagonal.*Global_H.*Omitted_Sigmae_x"))
