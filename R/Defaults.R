#' @export
DefaultArgsConfigOptim <- function(
  matParInit = NULL,
  argsPCMParamLowerLimit = NULL,
  argsPCMParamUpperLimit = NULL,
  numRunifInitVecParams = 1000,
  numGuessInitVecParams = 100,
  numJitterRootRegimeFit = 100, sdJitterRootRegimeFit = 0.5,
  numJitterAllRegimeFits = 100, sdJitterAllRegimeFits = 0.5,
  numCallsOptim = 10,
  control = list(fnscale = -1) ) {

  as.list(environment())
}

