# Copyright 2016-2020 Venelin Mitov
#
# This file is part of PCMFit.
#
# PCMFit is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PCMFit is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PCMFit.  If not, see <http://www.gnu.org/licenses/>.

#' Construct a MGPM object
#'
#' MGPM stays for Mixed Gaussian Phylogenetic Model. A MGPM object represents a
#' list of all variables that are subject to change during a probabilistic
#' inference of such a model for a given tree and trait data. By contrast, the
#' objects that change rarely (or not at all) during the inference, such as the
#' phylogetic tree, the trait data, the numerical limits for the different model
#' parameters, etc, are stored in a \code{\link{MGPMEnvironment}} object. By
#' convention, the tree in a MGPM is a \code{\link{PCMTree}} object
#' with node labels set to the character string representation of the corresponding
#' integer node ids. A \code{\link{PCMTree}} object extends the class
#' \code{phylo} from the R-package \code{ape} with additional fields specifying a
#' partition and a regime assignment (coloring) of the tree. The parts are named
#' as the character string labels of their corresponding root nodes. A regime,
#' i.e. a color, can be assigned to one or several non-neighboring parts and is
#' named as the corresponding part, which's name has the minimal integer value
#' among all parts, to which this regime is assigned.
#' Following this convention, a MGPM can be encoded as a numerical vector:
#' \deqn{\vec{s}=(K, n_2,...,n_{K+1}, l_2,...,l_{K+1}, r_2,...,r_{K+1}, m_1,...,m_R, v_1,...,v_P)^T}
#' The vector elements are described as follows:
#' \describe{
#' \item{\eqn{K}: }{number of shifts;}
#' \item{\eqn{(n_2,...,n_{K+1})^T}: }{integers representing labels of shift-nodes.
#' The shifts in the model occur at points within the branches leading to the
#' shift nodes in tip-ward direction. The corresponding locations of these
#' points are specified by \eqn{(l_2,...,l_{K+1})^T}.
#' The integers \eqn{(n_2,...,n_{K+1})^T} should be ordered in increasing order.}
#' \item{\eqn{(l_2,...,l_{K+1})^T}: }{offsets of the shift points measured as
#' fractions of the branch lengths in tip-ward direction. }
#' \item{\eqn{(r_2,...,r_{K+1})^T}: }{regime index vector. This is an integer
#' vector with elements among \eqn{(0,n_2,...,n_{K+1})^T}, indicating the
#' regime associated with each part in the tree. The regimes are named as the
#' shift nodes, except for 0, which corresponds to the root part (and regime).
#' The regime \eqn{r_1 = 0} is always present and, therefore,
#' omitted. It is possible to have lumped regimes, that is, different parts
#' of the tree having the same regime. This regime-lumping must obey the
#' following rules:
#' \enumerate{
#' \item neighbor parts cannot be lumped, i.e. they must differ by regime.
#' Two parts originating at nodes \eqn{n_i} and \eqn{n_j} in the tree are
#' called neighbor parts if they are separated solely by \eqn{n_i} or by \eqn{n_j};
#' \item to resolve the conflict between the shift nodes of the different parts
#' covered by a lumped regime, it is established that the name of a
#' regime assigned to several parts must equal the name of the part which is a
#' smaller integer number (remember that parts are named by the id's of their
#' root nodes, which, by convention are integer numbers).}
#' }
#' \item{\eqn{(m_1,...,m_R)^T}: }{model type assignment to the unique regimes.
#' This is an integer vector with elements between 1 and M, M denoting the
#' total number of model types possible. Each element corresponds to an
#' element in \code{sort(unique(c(0,r_2,...,r_{K+1})))}}
#' \item{\eqn{(v_1,...,v_P)^T}: }{real numbers passed to
#' \code{\link{PCMParamLoadOrStore}}. This is a vectorized form of the model
#' parameters.}
#' }
#' This function constructs a MGPM object either for the phylo object
#' \code{get("tree", env)} or for a subtree of this object defined by the nested
#' ED-expression \code{treeEDExpr}.
#'
#' @param env a \code{\link{MGPMEnvironment}} object.
#' @param K a single non-negative integer number denoting the number of shifts
#' in the model. Alternatively, this can be a double vector equal to the
#' concatenation of the double-converted arguments \code{K, n, l, r, m, v}. In
#' this case, the arguments \code{n, l, r, m, v} will be ignored. To check that
#' this is the case, the length of the argument K is checked for being bigger
#' than 1. This makes sense because the shortest possible concatenation of
#' \code{K, n, l, r, m, v} is of length 2. Setting this parameter to \code{NULL}
#' will cause it to be set at random with an integer in the range
#' \code{[0, M-N-1L]}, where \code{M-N-1L} is the number of internal nodes in
#' the tree. Default: \code{0L}.
#' @param n an integer vector of length \code{K} containing (integer) labels
#' of internal and/or tip nodes in \code{get("tree", env)} or the subtree evaluated by
#' \code{\link{PCMTreeEvalNestedEDxOnTree}(treeEDExpr, get("tree", env))}.
#' This argument is ignored if \code{length(K) > 1}. By default, or when
#' set to \code{NULL}, this argument will be set to
#' a random sample of \code{K} internal nodes of the tree obtained after
#' evaluating \code{treeEDExpr} on \code{get("tree", env)}. See also the argument
#' \code{treeEDExpr}.
#' @param l a double vector of length \code{K} with elements in [0,1]. Ignored if
#' \code{length(K) > 1}. When set to \code{NULL}, this argument will be drawn
#' from a uniform distribution. Default value: \code{rep(0.0, K[[1L]])}.
#' @param r an integer vector of length \code{(K)} containing regime names,
#' satisfying the rules in the description. If set to \code{NULL}, this argument is
#' assigned at random by respecting the rules in the description. Ignored if
#' \code{length(K) > 1}. Default value: \code{n}.
#' @param m an integer vector of length R, where R is the number of unique
#' regimes, specifying model type mapping for each regime. If set to \code{NULL},
#' this argument will be assigned at random. Ignored if
#' \code{length(K) > 1}. Default: \code{rep(1L, 1L + length(unique(r)))}.
#' @param v a double vector of length P where P is the number of variable
#' parameters of the MixedGaussian model corresponding to this MGPM object. If
#' set to \code{NULL}, this argument will be assigned at random using
#' the code \code{with(env, PCMParamRandomVecParams(model))}, where \code{model}
#' is the created \code{\link{MixedGaussian}} object corresponding to
#' \code{n, l, r, m, env} and \code{treeEDExpr}.
#' This argument will be ignored if \code{length(K) > 1}. By default this is set
#' to \code{NA_real_}, which results in a conversion to a vector of 0's.
#' @param treeEDExpr a character string (default "tree") corresponding to
#' a nested ED-expression for generating the MGPM object on a part of
#' \code{get("tree", env)}. See \code{\link{PCMTreeEvalNestedEDxOnTree}}.
#'
#' @return an object of S3 class 'MGPM'.
#' @examples
#' # The PCMBase package comes with a collection of simulated objects, which we
#' # can use for the examples.
#' library(PCMBase)
#'
#' tree <- PCMBaseTestObjects$tree.ab
#' X <- PCMBaseTestObjects$traits.ab.123[1:2, ]
#'
#' # Create a MGPM environment
#' env <- MGPMEnvironment(
#'   X, tree,
#'   model = MixedGaussian(
#'     k = 2,
#'     modelTypes = MGPMDefaultModelTypes(),
#'     mapping = structure(1:6, names = LETTERS[1:6]),
#'     Sigmae_x = Args_MixedGaussian_MGPMDefaultModelTypes()$Sigmae_x))
#'
#' mgpmDefault <- MGPM(env)
#'
#' stopifnot(MGPMPosK(mgpmDefault) == 1L)
#' stopifnot(identical(MGPMPosn(mgpmDefault), integer(0)))
#' stopifnot(identical(MGPMPosl(mgpmDefault), integer(0)))
#' stopifnot(identical(MGPMPosr(mgpmDefault), integer(0)))
#' stopifnot(identical(MGPMPosm(mgpmDefault), 2L))
#' stopifnot(identical(MGPMPosv(mgpmDefault), as.integer(3:6)))
#'
#' # A random MGPM
#' mgpmRandom <- MGPM(env, K = 6, n = NULL, l = rep(0, 6), r = NULL, m = NULL, v = NULL)
#' \donttest{
#' PCMTreePlot(
#'   mgpmRandom$tree) +
#'   ggtree::geom_nodelab(angle = 45) + ggtree::geom_tiplab(angle = 45)
#' }
#'
#'
#' mgpm <- MGPM(env,
#'   K = 4,                    # number of shifts
#'   n = c(8, 46, 52, 73),     # n_2, ..., n_{K+1}: shift nodes
#'   l = c(0, 0.4, 0.1, 0.05), # l_2, ..., l_{K+1}: offsets of the shift points
#'                             # relative to the beginnings of shift branches in
#'                             # tip-ward direction.
#'   r = c(0, 46, 52, 0),      # r_2,...,r_{K+1}: regime indices corresponding
#'                             # to the shifts. The regime of the part starting
#'                             # at the root node (41) is, by convention,
#'                             # always named 0 and is not included. The regime
#'                             # for the parts 8 and 73 is again 0, meaning that
#'                             # these regimes are lumped with the root regime.
#'                             # So the number of regimes are R=3, although there
#'                             # are 5 different parts in the tree.
#'   m = c(1, 5, 2)            # model type mapping for the R regimes.
#'                             # Note that these correspond to
#'                             # sort(unique(c(0, r_2, ..., r_{K+1}))).
#'                             # Hence, model type 1 corresponds to the
#'                             # single-branch regime (ending at tip 2), model
#'                             # type 5 correspods to the root regime (41),
#'                             # model type 2 corresponds to 46 and model type 3
#'                             # corresponds to 52.
#'   )
#'
#' \donttest{
#' PCMTreePlot(
#'   mgpm$tree) +
#'   ggtree::geom_nodelab(angle = 45) + ggtree::geom_tiplab(angle = 45)
#' }
#'
#' mgpm$n
#'
#' # the names correspond to the shift nodes.
#' mgpm$l
#' mgpm$r
#'
#' # The names of mgpm$m denote the regimes, while the values denote the model
#' # types:
#' mgpm$m
#' PCMTreeGetPartRegimes(mgpm$tree)
#'
#' # This should fail because the shift nodes are not sorted in increasing order.
#' error <- try(
#'   MGPM(env, 4,
#'        c(8, 52, 46, 73),
#'        c(0, 0.1, 0.4, 0.0),
#'        c(0, 52, 46, 0),
#'        c(1, 5, 2)), silent = TRUE)
#' stopifnot(inherits(error, "try-error"))
#'
#' # At this stage, mgpm's parameters are all equal to 0.
#' # Don't do this, because it does not execute PCMParamRandomVecParams in the environment env.
#' mgpmVec <- as.vector(mgpm)
#' mgpmVec[MGPMPosv(mgpm)] <- PCMParamRandomVecParams(mgpm$model)
#'
#' # Do this way:
#' mgpm2 <- MGPM(attr(mgpm, "env"), mgpm$K, mgpm$n, mgpm$l, mgpm$r, mgpm$m, NULL)
#'
#' stopifnot(identical(mgpm$n, mgpm2$n))
#' stopifnot(identical(mgpm$l, mgpm2$l))
#' stopifnot(identical(mgpm$r, mgpm2$r))
#' stopifnot(identical(mgpm$m, mgpm2$m))
#'
#' # These two should be different values but same length:
#' cbind(mgpm$v, mgpm2$v)
#'
#' mgpm3 <- MGPM(env, c(4, 8, 46, 52, 73, 0, 0.4, 0.1, 0.05, 0, 46, 52, 0, 1, 5, 2))
#' stopifnot(identical(mgpm$n, mgpm3$n))
#' stopifnot(identical(mgpm$l, mgpm3$l))
#' stopifnot(identical(mgpm$r, mgpm3$r))
#' stopifnot(identical(mgpm$m, mgpm3$m))
#'
#' mgpm4 <- MGPM(env,
#'   K = 1,
#'   n = c(62),
#'   l = c(0.2),
#'   r = c(62),
#'   m = c(1, 5),
#'   treeEDExpr = "D(E(tree,52),73)")
#' PCMTreeGetPartRegimes(mgpm4$tree)
#'
#' \donttest{
#' PCMTreePlot(
#'   mgpm4$tree) +
#'   ggtree::geom_nodelab(angle = 45) + ggtree::geom_tiplab(angle = 45)
#' }
#' @seealso \code{\link{MGPMEnvironment}} \code{\link{MGPMPosK}} \code{\link{PCMBase::MixedGaussian}}
#'
#' @import PCMBase
#' @export
MGPM <- function(
  env,
  K = 0L,
  n = as.integer(sample(PCMTreeGetNodeLabels(PCMTreeEvalNestedEDxOnTree(treeEDExpr, get("tree", env))), K)),
  l = rep(0.0, K[[1L]]),
  r = n,
  m = rep(1L, 1L + length(unique(r))),
  v = NA_real_,
  treeEDExpr = "tree") {

  if(!is.MGPMEnvironment(env)) {
    stop("MGPM:: env should be an MGPMEnvironment object.")
  }

  ##................................................................
  ##
  ##  Create the tree for the MGPM object
  ##
  ##................................................................

  # The node labels in get("tree", env) must be 1,2,...,M, with N+1 being the root and
  # 1,...,N being the tips.
  tree <- get("tree", env)

  if(!identical(PCMTreeGetLabels(tree),
                as.character(seq_len(PCMTreeNumNodes(tree))))) {
    stop(paste0(
      'MGPM:: The node labels in get("tree", env) must be 1,2,...,M, with N+1 being the root ',
      'and 1,...,N being the tips.'))
  }

  if(!is.null(treeEDExpr) && !identical(treeEDExpr, "tree")) {
    tree <- try(PCMTreeEvalNestedEDxOnTree(treeEDExpr, tree), silent = TRUE)
    if(inherits(tree, "try-error")) {
      stop(paste0("MGPM:: Error while processing treeEDExpr:", tree))
    }
  }

  # number of tips in the tree
  N <- PCMTreeNumTips(tree)
  M <- PCMTreeNumNodes(tree)


  ##................................................................
  ##
  ##  K: number of shifts or a vector containing K,n,l,r,m,v.
  ##
  ##................................................................

  if(is.null(K)) {
    # Generate a random K
    K <- sample(0:(M-N-1L), 1L)
  }

  if(length(K) > 1L) {
    # s: a numerical vector containing K,n,l,r,m,v
    s <- as.double(K)
  } else {
    # Arguments K,n,l,r,m,v are specified or are set to NULL.
    s <- NULL
  }

  if(is.null(s)) {
    # do nothing
  } else {
    # the first element is the number of shifts
    K <- as.integer(s[1L])
  }

  # avoid "no visible binding NOTE":
  i <- r2 <- NULL

  ##................................................................
  ##
  ##  n: shift nodes;
  ##  l: shift offsets from the beginning of branches (fractions)
  ##
  ##................................................................

  if(is.null(s)) {
    if(is.null(n)) {
      # Generate random shift nodes
      n <- sort(as.integer(sample(PCMTreeGetNodeLabels(tree), K)))
    }
    if(is.null(l)) {
      # Generate random shift offsets
      l <- runif(K)
    }
    if(length(n) != K || length(l) != K) {
      stop(
        paste0("MGPM:: Some of (n, l) is/are of length different than K(", K, ")."))
    }
  } else {
    n <- as.integer(s[1L + seq_len(K)])
    l <-  as.double(s[1L + K + seq_len(K)])
  }

  # n, l and r should be ordered in the increasing order of n
  ordern <- order(n)
  if(!identical(ordern, seq_len(K))) {
    stop("MGPM:: the shift nodes (n) should be ordered in increasing order.")
  }

  ##................................................................
  ##
  ##  Add singleton nodes to the tree where needed.
  ##
  ##................................................................

  # Insert singleton nodes in branches where the shift-point is different from the
  # beginning of the branch.
  dtShifts <- data.table(
    n = c(0L, as.integer(n)),
    l = as.double(c(0.0, l)))

  if(isTRUE(any(dtShifts[-1L, l>1 | l<0]))) {
    stop("MGPM:: elements of l should be in the interval [0,1].")
  }

  minBranchLength <- mget("minBranchLength", env, ifnotfound = .05)[[1L]]

  # named branch lengths
  blen <- structure(tree$edge.length, names = as.character(tree$edge[, 2L]))

  dtShifts[, nId:=N+1L]
  dtShifts[-1L, nId:=PCMTreeMatchLabels(tree, as.character(n))]

  dtShifts[-1L, len:=blen[as.character(nId)]]
  # absolute offset in root-ward direction (needed for PCMTreeInsertSingletons)
  dtShifts[-1L, lPos:=len-l*len]
  dtShifts[, needASingleton:=FALSE]
  dtShifts[-1L, needASingleton:=(lPos >= minBranchLength & len-lPos >= minBranchLength)]

  if(sum(dtShifts[, needASingleton]) > 0L) {
    tree <- PCMTreeInsertSingletons(
      tree,
      dtShifts[needASingleton == TRUE, nId],
      dtShifts[needASingleton == TRUE, lPos])

    # update the id's of the shift nodes
    dtShifts[1L, nId:=PCMTreeNumTips(tree) + 1L]
    dtShifts[-1L, nId:=PCMTreeMatchLabels(tree, as.character(n))]
  }

  PCMTreeSetPartition(tree, dtShifts[, nId])
  dtShifts[, part:=PCMTreeGetPartsForNodes(tree, nodes = nId)]

  ##................................................................
  ##
  ##  r: regime for the parts
  ##
  ##  Currently each part in tree has its own regime, which is named in the same
  ##  way as the part.
  ##  We now process the r argument which specifies which parts will be "lumped",
  ##  i.e. share the same regime.
  ##  The following rules must be respected:
  ##  1. Neighbor parts don't have the same regime
  ##  2. Regimes that cover different parts must be named as the name of the part,
  ##  which is the smaller integer.
  ##
  ##................................................................


  dtShifts[-1L, nParentId:=sapply(nId, PCMTreeGetParent, tree = tree)]
  dtShifts[-1L, partParent:=PCMTreeGetPartsForNodes(tree, nodes = nParentId)]
  dtShifts[, daughterParts:=lapply(.I, function(i) part[partParent == part[i]][-1L])]

  # make search over n binary
  setkey(dtShifts, n)

  rIsRandom <- FALSE
  if(is.null(s)) {
    if(is.null(r)) {
      r <- n
      # We generate r at random while respecting the rules
      rIsRandom <- TRUE
    } else {
      r <- try(as.integer(r), silent = TRUE)
      if(inherits(r, "try-error") || length(r) != K) {
        stop("MGPM:: r should be an integer vector of length K.")
      }
    }
  } else {
    r <- try(as.integer(s[1L + K + K + seq_len(K)]), silent = TRUE)
    if(inherits(r, "try-error")) {
      stop(paste0("MGPM:: error while parsing the regime vector r: ", r))
    }
  }

  rootRegime = structure(0L, names = PCMTreeGetRootLabel(tree))
  names(r) <- as.character(n)

  for(i in seq_len(K)) {
    # Candidate regimes for shift i: By rule 2, the candidates are either i, 0,
    # or any of the regimes
    # of shift nodes that precede n[i], such that the parent of n[i] does not
    # belong to their part.
    shiftNode <- n[i]
    regimePartParent <- c(rootRegime, r)[dtShifts[n==shiftNode, partParent]]
    regimesDaughterParts <- c(rootRegime, r)[dtShifts[n==shiftNode,
                                                      daughterParts[[1L]]]]

    candRegimes <- unname(setdiff(
      c(0L, r[seq_len(i)]), union(regimePartParent, regimesDaughterParts)))

    if(rIsRandom) {
      if(length(candRegimes) > 1L) {
        r[i] <- sample(candRegimes, 1L)
      } else {
        r[i] <- candRegimes[1L]
      }
    } else {
      if( !(r[i] %in% candRegimes) ) {
        stop(paste0(
          "MGPM:: The specified vector r violates the rule that neighbor parts ",
          "should not be lumped."))
      }
    }
  }

  dtShifts[, r:=c(0L, r)]

  # Check that lumped regimes are correctly specified and set their names to the
  # smallest shift-node
  dtShifts[, r2:=min(n), by=r]
  if(!dtShifts[, identical(r, r2)]) {
    stop(paste0(
      "MGPM:: Except for the root-part regime, which is to be named 0, each ",
      "regime lumping several parts of the tree should be named ",
      "as the smallest id of a shift-node among the shift-nodes that are at ",
      "the root of the lumped parts."))
  }


  if(length(dtShifts[, unique(c(n, r))]) != K+1L) {
    stop(paste0(
      "MGPM:: There are duplicated root or shift-nodes or some of the ",
      "regime names may not be among the root and the shift-nodes. \n",
      "root node and shift nodes: ", toString(dtShifts[, n]), "\n",
      "regimes: ", toString(dtShifts[, r])))
  }

  PCMTreeSetPartRegimes(tree, dtShifts[, structure(as.character(r), names = part)])

  ##................................................................
  ##
  ##  m: model types mapped to the regimes
  ##
  ##................................................................

  # unique regime names
  regimeNames <- dtShifts[, unique(r)]

  # Number of unique regimes
  R <- length(regimeNames)

  if(is.null(s)) {
    if(is.null(m)) {
      m <- sample(seq_along(get("modelTypeNames", env)), R, replace = TRUE)
    }
    m <- as.integer(m)
    if(length(m) != R) {
      stop(paste0(
        "MGPM:: the length of m (", length(m),
        ") should equal the number of unique regimes (", R, ")."))
    }
  } else {
    # model types associated with the regimes
    if(length(s) >= 1L + K + K + K + R) {
      m <- as.integer(s[1L + K + K + K + seq_len(R)])
    } else {
      stop("MGPM:: not enough model types supplied in input (s).")
    }
  }
  names(m) <- regimeNames
  env$mNewModel <- m

  ##................................................................
  ##
  ##  model: The MixedGaussian model
  ##
  ##................................................................

  with(env, newModel <- do.call(
    MixedGaussian,
    c(list(k = k,
           modelTypes = modelTemplate[modelTypeNames],
           mapping = mNewModel),
      attr(modelTemplate, "spec")[setdiff(names(attr(modelTemplate, "spec")), modelTypeNames)])))

  # model <- do.call(
  #   MixedGaussian,
  #   c(list(k = get("k", env),
  #          modelTypes = get("modelTemplate", env)[get("modelTypeNames", env)],
  #          mapping = m),
  #     attr(get("modelTemplate", env), "spec")[
  #       setdiff(names(attr(get("modelTemplate", env), "spec")), get("modelTypeNames", env))]))

  # Load the parameter values into the model (the number of parameters is
  # stored in P)
  P <- as.integer(attr(env$newModel, "p"))

  ##................................................................
  ##
  ##  v: The vector of model parameter values
  ##
  ##................................................................

  if(is.null(s)) {
    if(is.null(v)) {
      # Generate a random v according to the parameter limits. This is done by evaluating
      # PCMParamRandomVecParams(model) in env, in order to use any specific settings such
      # as S3 methods for the functions PCMBase::PCMParamLowerLimit and
      # PCMBase::PCMParamUpperLimit. However.
      # The object model must exist in the environme
      v <- with(env, PCMParamRandomVecParams(newModel))
    } else if(is.na(v)) {
      v <- double(P)
    } else if(length(v) != P) {
      stop(paste0(
        "MGPM:: if v is specified, the length of v (", length(v),
        ") should equal the number of variable parameters of the model (", P, ")"))
    } else {
      # Do nothing, i.e. use the specified v
    }
  } else {
    v <- double(P)
    if(length(s) >= 1L + K + K + K + R + P) {
      v <- s[-seq_len(1L + K + K + K + R)]
    }
  }

  model <- env$newModel
  P1 <- PCMParamLoadOrStore(model, v, offset = 0, k = get("k", env), R = R, load = TRUE)
  if(P != P1) {
    stop(paste0(
      "MGPM:: something is wrong with the attribute p of the model ",
      "and the number of parameters returned by PCMParamLoadOrStore. ",
      "This could be a bug."))
  }

  # Clean newModel and mNewModel from env
  rm("newModel", "mNewModel", pos = env)

  ##................................................................
  ##
  ##  MGPM object
  ##
  ##................................................................

  res <- structure(list(
    K = K, R = R, P = P,
    n = dtShifts[-1L, n],
    l = structure(dtShifts[-1L, l], names = dtShifts[-1L, n]),
    r = structure(dtShifts[-1L, r], names = dtShifts[-1L, n]),
    m = m,
    v = v,
    model = model,
    tree = tree),
    class = "MGPM")
  attr(res, "env") <- env
  attr(res, "treeEDExpr") <- treeEDExpr

  res
}


#' Check if an object is of S3 class 'MGPM'
#' @param o an object.
#' @return a logical.
#' @export
is.MGPM <- function(o) {
  inherits(o, "MGPM")
}

#' If necessary, convert an object to an MGPM object.
#' This is an S3 generic function.
#' @param o an object.
#' @param env a \code{\link{MGPMEnvironment}} object. Default: \code{attr(o, "env")}.
#' @param treeEDExpr a character string denoting a ED-expression.
#' Default: \code{attr(o, "treeEDExpr")}.
#' \code{\link{PCMBase::PCMTreeEvalNestedEDxOnTree()}}.
#'
#' @return an MGPM object corresponding to \code{o}. If \code{o} is already
#' a \code{\link{MGPM}} object, it is returned as is. Otherwise, a
#' \code{\link{MGPM}} object is constructed from \code{o} and,
#' optionally \code{env} and \code{treeEDExpr}.
#' @seealso \code{\link{MGPM}}
#' @export
as.MGPM <- function(o, env = attr(o, "env"), treeEDExpr = attr(o, "treeEDExpr")) {
  UseMethod("as.MGPM", o)
}

#' @export
as.MGPM.default <- function(o, env = attr(o, "env"), treeEDExpr = attr(o, "treeEDExpr")) {
  stop(paste0("as.MGPM.default:: no method defined for class ",
              toString(class(o)), "." ))
}

#' @export
as.MGPM.MGPM <- function(o, env = attr(o, "env"), treeEDExpr = attr(o, "treeEDExpr")) {
  o
}

#' @export
as.MGPM.double <- function(o, env = attr(o, "env"), treeEDExpr = attr(o, "treeEDExpr")) {
  MGPM(env, K = o, treeEDExpr = treeEDExpr)
}

#' @export
as.vector.MGPM <- function(x, mode = "any") {
  c(as.double(unname(x$K)),
    as.double(unname(x$n)),
    as.double(unname(x$l)),
    as.double(unname(x$r)),
    as.double(unname(x$m)),
    as.double(unname(x$v)))
}


#' @title Indices of different parts of an MGPM in an MGPMVector
#'
#' @description The different parts of an MGPM are described in
#' \code{\link{MGPMVector}}. \code{posK} returns the osition of the
#'  number of shifts (this is always equal to 1). See the section functions
#'  for the others.
#'
#' @param s a \code{\link{MGPM}} or a \code{\link{MGPMVector}}
#' object.
#' @param env a \code{\link{MGPMEnvironment}} object. This is needed only if
#'  \code{s} is \code{\link{MGPMVector}} object.
#'
#' @return an integer vector.
#'
#' @seealso \code{\link{MGPM}}
#' @export
MGPMPosK <- function(s, env = NULL) {
  1L
}

#' @describeIn MGPMPosK
#'
#' Positions of the shift nodes;
#'
#' @export
MGPMPosn <- function(s, env = NULL) {
  s <- as.MGPM(s, env)
  1L + seq_len(s$K)
}

#' @describeIn MGPMPosK
#'
#' Positions of the offsets of the shift-points in tip-ward direction
#' relative to the beginnings of the branches leading to shift nodes;
#'
#' @export
MGPMPosl <- function(s, env = NULL) {
  s <- as.MGPM(s, env)
  1L + s$K + seq_len(s$K)
}

#' @describeIn MGPMPosK
#'
#' Position of the regime id's for each shift node;
#'
#' @export
MGPMPosr <- function(s, env = NULL) {
  s <- as.MGPM(s, env)
  1L + s$K + s$K + seq_len(s$K)
}

#' @describeIn MGPMPosK
#'
#' Positions of the model type id's for each regime;
#'
#' @export
MGPMPosm <- function(s, env = NULL) {
  s <- as.MGPM(s, env)
  1L + s$K + s$K + s$K + seq_len(s$R)
}

#' @describeIn MGPMPosK
#'
#' Positions of the model parameters.
#'
#' @export
MGPMPosv <- function(s, env = NULL) {
  s <- as.MGPM(s, env)
  1L + s$K + s$K + s$K + s$R + seq_len(s$P)
}

