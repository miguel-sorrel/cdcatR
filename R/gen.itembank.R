#' Item bank generation
#'
#' This function can be used to generate an item bank.
#' The user can provide a Q-matrix or create one defining a set of arguments.
#' Item quality is sampled from a uniform distribution with mean = \emph{mean.IQ} and range = \emph{range.IQ}.
#' Item parameters are generated so that the monotonicity constraint is satisfied.
#'
#' @param Q Numeric matrix of length \emph{J} number of items x \emph{K} number of atributes. Q-matrix
#' @param gen.Q A list of arguments to generate a Q-matrix if \code{Q} is not provided. \code{J} = number of items (scalar numeric). \code{K} = number of attributes (scalar numeric). \code{propK.J} = numeric vector summing up to 1 that determines the proportion of 1-attribute, 2-attribute, ..., items. The length of \code{propK.J} determines the maximum number of attributes considered for an item (see \code{Examples} below). \code{nI} = Scalar numeric that sets the minimum number of identity matrices to be included in the Q-matrix. \code{minJ.K} = numeric vector of length \emph{K} that sets the minimum number of items measuring each attribute. \code{max.Kcor} = scalar numeric that sets the maximum positive correlation allowed between two attributes
#' @param mean.IQ Item discrimination (mean for the uniform distribution). \emph{mean.IQ} = \emph{P}(\strong{1}) - \emph{P}(\strong{0}) (Sorrel et al., 2017; Najera et al., in press). Must be a scalar numeric between 0 and 1
#' @param range.IQ Item discrimination (range for the uniform distribution). Must be a scalar numeric between 0 and 1
#' @param model A character vector of length \emph{J} with one model for each item, or a single value to be used for all items. The possible options include \code{"DINA"}, \code{"DINO"}, \code{"ACDM"}, and \code{"GDINA"}. One-attribute items will be coded in the output as \code{"GDINA"}
#' @param min.param Scalar numeric. Minimum value for the delta parameter of the principal effects of each attribute. Only usable if \code{model} = \code{"ACDM"} or \code{model} = \code{"GDINA"}
#' @param seed Scalar numeric. A scalar to use with \code{set.seed}
#'
#' @return \code{gen.itembank} returns an object of class \code{gen.itembank}.
#' \describe{
#' \item{simQ}{Generated Q-matrix (only if \code{gen.Q} arguments have been used)}
#' \item{simcatprob.parm}{A list of success probabilities for each latent group in each item}
#' \item{simdelta.parm}{A list of delta parameters for each item}
#' \item{check}{A list that contains the mean.IQ and range.IQ for the item bank so that users can check whether these values match the expected results}
#' \item{specifications}{A list that contains all the specifications}
#' }
#'
#' @references
#'
#' Najera, P., Sorrel, M. A., de la Torre, J., & Abad, F. J. (in press). Improving robustness in Q-matrix validation using an iterative and dynamic procedure. \emph{Applied Psychological Measurement}.
#'
#' Sorrel, M. A., Abad, F. J., Olea, J., de la Torre, J., & Barrada, J. R. (2017). Inferential item-fit evaluation in cognitive diagnosis modeling. \emph{Applied Psychological Measurement, 41}, 614-631.
#'
#' @examples
#'####################################
#'# Example 1.                       #
#'# Generate item bank providing a   #
#'# Q-matrix using the G-DINA model  #
#'####################################
#'
#' Q <- sim180GDINA$simQ
#' bank <- gen.itembank(Q = Q, mean.IQ = .70, range.IQ = .20, model = "GDINA")
#'
#'####################################
#'# Example 2.                       #
#'# Generate item bank providing a   #
#'# Q-matrix using multiple models   #
#'####################################
#'
#' Q <- sim180GDINA$simQ
#' K <- ncol(Q)
#' model <- sample(c("DINA", "DINO", "ACDM"), size = nrow(Q), replace = TRUE)
#' bank <- gen.itembank(Q = Q, mean.IQ = .70, range.IQ = .20, model = model)
#'
#'####################################
#'# Example 3.                       #
#'# Generate item bank without       #
#'# providing a Q-matrix (using      #
#'# gen.Q arguments)                 #
#'####################################
#'
#' bank <- gen.itembank(gen.Q = list(J = 150, K = 5, propK.J = c(0.4, 0.3, 0.2, 0.1),
#'                      nI = 3, minJ.K = 30, max.Kcor = 1),
#'                      mean.IQ = .80, range.IQ = .10, min.param = 0.1)
#'
#' @export
#'
gen.itembank <- function(Q = NULL,
                         gen.Q = list(J = NULL, K = NULL, propK.J = NULL, nI = 1, minJ.K = NULL, max.Kcor = 1),
                         mean.IQ, range.IQ, model = "GDINA", min.param = 0,
                         seed = NULL)
  {

  #----------------------------
  # Warning end error messages
  #----------------------------

  if(is.null(Q) & any(sapply(gen.Q, is.null))){stop("Q argument or all gen.Q arguments required.")}
  if(!is.null(Q) & all(!sapply(gen.Q, is.null))){warning("If Q argument is specified, gen.Q arguments are ignored.")}
  if(is.null(mean.IQ)){stop("mean.IQ argument missing, with no default.")}
  if(mean.IQ > 1 | mean.IQ < 0){stop("mean.IQ must be a value between 0 and 1.")}
  if(is.null(range.IQ)){stop("range.IQ argument missing, with no default.")}
  if(range.IQ > 1 | range.IQ < 0){stop("range.IQ must be a value between 0 and 1.")}
  if(is.null(Q)){
    if(is.null(gen.Q$J)){stop("gen.Q$J argument required.")}
    if(is.null(gen.Q$K)){stop("gen.Q$K argument required.")}
    if(is.null(gen.Q$propK.J)){stop("gen.Q$propK.J argument required.")}
    if(is.null(gen.Q$nI)){gen.Q$nI <- 1}
    if(is.null(gen.Q$minJ.K)){stop("gen.Q$minJ.K argument required.")}
    if(is.null(gen.Q$max.Kcor)){gen.Q$max.Kcor <- 1}
  }

  #----------------
  # Gather objects
  #----------------

  arg.range.IQ <- range.IQ
  range.IQ <- range.IQ/2

  if(!is.null(seed)){set.seed(seed)}
  if(is.null(seed)){seed <- sample(1:1000000, size = 1)}

  if(!is.null(Q)){
    sim.Q <- NULL
    arg.Q <- Q
    J <- nrow(Q)
    K <- ncol(Q)
    arg.model <- model
    if(length(model) == 1){model <- rep(model, J)}
  } else {
    arg.Q <- NULL
    J <- gen.Q$J
    K <- gen.Q$K
    sim.Q <- Q <- genQ(J, K, gen.Q$nI, gen.Q$propK.J, T, gen.Q$minJ.K, gen.Q$max.Kcor, seed)$Q
    arg.model <- model
    if(length(model) == 1){model <- rep(model, J)}
  }
  min.P0 <- round(((1 - mean.IQ) / 2) - (range.IQ / 2), 10)
  max.P0 <- round(((1 - mean.IQ) / 2) + (range.IQ / 2), 10)
  min.P1 <- round(((1 - ((1 - mean.IQ) / 2))) - (range.IQ / 2), 10)
  max.P1 <- round(((1 - ((1 - mean.IQ) / 2))) + (range.IQ / 2), 10)
  if(min.P0 < 0){
    warning("mean.IQ and range.IQ may give probabilities below 0 for P0. Probabilities coherced to be 0 at least.")
    min.P0 <- 0
  }
  if(max.P1 > 1){
    warning("mean.IQ and range.IQ may give probabilities above 1 for P1. Probabilities coherced to be 1 at most.")
    max.P1 <- 1
  }
  if(max.P0 > min.P1){
    warning("mean.IQ and range.IQ may give to P0 higher probabilities than to P1. Probabilities coherced to be always higher for P1.")
    diff <- max.P0 - min.P1
    max.P0 <- max.P0 - diff / 2
    min.P1 <- min.P1 + diff / 2
  }
  P <- cbind(P0 = runif(J, min.P0, max.P0), P1 = runif(J, min.P1, max.P1))

  #--------------------------
  # Generate item parameters
  #--------------------------
  catprob.parm <- delta.parm <- list()
  for(j in 1:J){
    k.j <- sum(Q[j,])
    lg.j <- 2^sum(Q[j,])
    param <- switch(model[j],
                    DINA = list(catprob.parm = c(rep(P[j,"P0"], lg.j - 1), P[j,"P1"]),
                                delta.parm = c(P[j,"P0"], P[j,"P1"] - P[j,"P0"])),
                    DINO = list(catprob.parm = c(P[j,"P0"], rep(P[j,"P1"], lg.j - 1)),
                                delta.parm = c(P[j,"P0"], P[j,"P1"] - P[j,"P0"])),
                    ACDM = gen.ACDMparam(P[j,"P0"], P[j,"P1"], k.j, min.param),
                    GDINA = gen.GDINAparam(P[j,"P0"], P[j,"P1"], k.j, min.param))
    catprob.parm[[j]] <- param$catprob.parm
    delta.parm[[j]] <- param$delta.parm
    names(catprob.parm[[j]]) <- names(delta.parm[[j]]) <- NULL
  }

  check <- list()
  check[["mean.IQ.bank"]] <-
    mean(P[,2] - P[,1])
  check[["range.IQ.bank"]] <-
    max(P[,2]-P[,1])-min(P[,2]-P[,1])
  if(is.null(arg.Q)){
  check[["propK.J"]] <- table(rowSums(sim.Q))/nrow(sim.Q)
  tmp <- cor(sim.Q)
  diag(tmp) <- NA
  check[["max.Kcor"]] <- max(tmp, na.rm = TRUE)}

  #----------------
  # Return results
  #----------------

  res <- list()
  res$simQ <- sim.Q
  res$simcatprob.parm <- catprob.parm
  res$simdelta.parm <- delta.parm
  res$check <- check
  model[which(rowSums(Q)==1)] <- "GDINA"
  res$specifications <- list(Q = arg.Q, gen.Q = gen.Q, mean.IQ = mean.IQ, range.IQ = arg.range.IQ,
                             model = model, min.param = min.param, seed = seed)

  class(res) <- "gen.itembank"
  return(res)
}
