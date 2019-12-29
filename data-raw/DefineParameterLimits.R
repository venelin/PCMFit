# lower limits for models A and B
PCMParamLowerLimit.BM <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  if(is.Global(o$Sigma_x)) {
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- -.0
    }
  } else {
    if(!is.Diagonal(o$Sigma_x)) {
      for(r in seq_len(R)) {
        o$Sigma_x[1, 2, r] <- -.0
      }
    }
  }
  o
}

# upper limits for models A and B
PCMParamUpperLimit.BM <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[1, 1] <- o$Sigma_x[2, 2] <- 1.0
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- 1.0
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[1, 1, r] <- o$Sigma_x[2, 2, r] <- 1.0
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2, r] <- 1.0
      }
    }
  }
  o
}

# lower limits for models C, ..., F.
PCMParamLowerLimit.OU <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  if(is.Global(o$Theta)) {
    o$Theta[1] <- 0.0
    o$Theta[2] <- -1.2
  } else {
    for(r in seq_len(R)) {
      o$Theta[1, r] <- 0.0
      o$Theta[2, r] <- -1.2
    }
  }
  if(is.Global(o$Sigma_x)) {
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- -.0
    }
  } else {
    if(!is.Diagonal(o$Sigma_x)) {
      for(r in seq_len(R)) {
        o$Sigma_x[1, 2, r] <- -.0
      }
    }
  }
  o
}

# upper limits for models C, ..., F.
PCMParamUpperLimit.OU <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  if(is.Global(o$Theta)) {
    o$Theta[1] <- 7.8
    o$Theta[2] <- 4.2
  } else {
    for(r in seq_len(R)) {
      o$Theta[1, r] <- 7.8
      o$Theta[2, r] <- 4.2
    }
  }
  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[1, 1] <- o$Sigma_x[2, 2] <- 1.0
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- 1.0
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[1, 1, r] <- o$Sigma_x[2, 2, r] <- 1.0
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2, r] <- 1.0
      }
    }
  }
  o
}
