#' Item bank generation
#'
#' This function can be used to generate a calibration and a validation sample of responses to an item bank.
#' The user can provide a Q-matrix or create one defining a set of arguments.
#' Item quality is sampled from a uniform distribution with mean = \emph{IQ} and variance = \emph{VAR}.
#' Item parameters are generated with monotonicity constraint.
#' Data are generated following either a uniform, multivariate normal, or higher-order distribution.
#'
#' @param N.c Sample size for the calibration sample
#' @param N.v Sample size for the validation sample. If NULL (by default), no validation data are provided
#' @param Q Q-matrix
#' @param gen.Q A list of arguments to generate a Q-matrix if \code{Q} is not provided. \code{J}: number of items. \code{K}: number of attributes. \code{propK.J}: vector summing up to 1 that determines the proportion of 1-attribute, 2-attribute, ..., items (See \code{Examples}). \code{nI}: minimum number of identity matrices to be included in the Q-matrix. \code{minJ.K}: minimum number of items measuring each attribute. \code{max.Kcor}: maximum correlation allowed between two attributes
#' @param IQ Item discrimination (mean for the uniform distribution). \emph{IQ} = \emph{P}(\strong{1}) - \emph{P}(\strong{0}) (Sorrel, Abad, Olea, de la Torre, and Barrada, 2017; NÃ¡jera, Sorrel, de la Torre, and Abad, in press)
#' @param VAR Item discrimination (variance for the uniform distribution)
#' @param model A character vector with one model for each item, or a single value to be used for all items. The possible options include \code{"DINA"}, \code{"DINO"}, \code{"ACDM"}, and \code{"GDINA"}. One-attribute items will be coded as \code{"GDINA"}
#' @param min.param Minimum value for the delta parameter of the principal effects of each attribute. Only usable if \code{model} = \code{"ACDM"} or \code{model} = \code{"GDINA"}
#' @param gen.att A list of arguments for the attribute profiles generation. Uniform attribute distribution by default. See \code{simGDINA} function of package \code{GDINA}
#' @param seed A scalar to use with \code{set.seed}
#'
#' @return \code{gen.itembank} returns an object of class \code{gen.itembank}.
#' \describe{
#' \item{simdat.c}{A matrix with the calibration dataset}
#' \item{simalpha.c}{A matrix with the alpha patterns for the calibration dataset}
#' \item{simdat.v}{A matrix with the validation dataset (only if \code{N.v} is specified)}
#' \item{simalpha.v}{A matrix with the alpha patterns for the validation dataset (only if \code{N.v} is specified)}
#' \item{simQ}{Q-matrix}
#' \item{simcatprob.parm}{A list of non-zero categories success probabilities for each latent group in each item}
#' \item{simdelta.parm}{A list of delta parameters for each item}
#' \item{specifications}{A list that contains all the specifications}
#' }
#'
#' @references
#'
#' Kaplan, M., de la Torre, J., & Barrada, J. R. (2015). New item selection methods for cognitive diagnosis computerized adaptive testing. \emph{Applied Psychological Measurement, 39}, 167-188.
#'
#' Najera, P., Sorrel, M. A., de la Torre, J., & Abad, F. J. (in press). Improving robustness in Q-matrix validation using an iterative and dynamic procedure. \emph{Applied Psychological Measurement}.
#'
#' Sorrel, M. A., Abad, F. J., Olea, J., de la Torre, J., & Barrada, J. R. (2017). Inferential item-fit evaluation in cognitive diagnosis modeling. \emph{Applied Psychological Measurement, 41}, 614-631.
#'
#' @examples
#'####################################
#'# Example 1.                       #
#'# Generate calibration sample      #
#'# and item parameters providing    #
#'# a Q-matrix, with G-DINA model    #
#'# and uniform attribute            #
#'# distribution                     #
#'####################################
#'
#' Q <- sim180GDINA$simQ
#' model <- rep("GDINA", each = nrow(Q))
#' IQ <- .70 # P(1), IQ = Low item quality in Kaplan, de la Torre & Barrada (2015)
#' VAR <- .10 # High variance in Kaplan et al. (2015)
#' bank <- gen.itembank(N.c = 1000, Q = Q, IQ = IQ, VAR = VAR, model = model)
#'
#'####################################
#'# Example 2.                       #
#'# Generate calibration sample      #
#'# and item parameters providing    #
#'# a Q-matrix, with different       #
#'# models and higher-order          #
#'# attribute distribution           #
#'####################################
#'
#' Q <- sim180GDINA$simQ
#' K <- ncol(Q)
#' model <- sample(c("DINA", "DINO", "ACDM"), size = nrow(Q), replace = TRUE)
#' IQ <- .70
#' VAR <- .10
#' bank <- gen.itembank(N.c = 1000, Q = Q, IQ = IQ, VAR = VAR, model = model,
#'                      gen.att = list(att.dist = "higher.order",
#'                       higher.order.a = runif(K, 0.8, 1.3),
#'                       higher.order.b = seq(-2, 2, length.out = K)))
#'
#'####################################
#'# Example 3.                       #
#'# Generate calibration and         #
#'# validation samples and item      #
#'# parameters without providing a   #
#'# Q-matrix (using gen.Q arguments) #
#'####################################
#'
#' J <- 150
#' K <- 5
#' propK.J <- c(0.5, 0.3, 0, 0.2)
#' nI <- 3
#' minJ.K <- 50
#' max.Kcor <- 0.5
#' IQ <- .80
#' VAR <- .05
#' min.param <- 0.1
#' bank <- gen.itembank(N.c = 1000, N.v = 2000,
#'                      gen.Q = list(J = J, K = K, propK.J = propK.J,
#'                      nI = nI, minJ.K = minJ.K, max.Kcor = max.Kcor),
#'                      IQ = IQ, VAR = VAR, min.param = min.param)
#'
#' @export
#'
gen.itembank <- function(N.c, N.v = NULL, Q = NULL,
                         gen.Q = list(J = NULL, K = NULL, propK.J = NULL, nI = NULL, minJ.K = NULL, max.Kcor = NULL),
                         IQ, VAR, model = "GDINA", min.param = 0,
                         gen.att = list(att.dist = "uniform", mvnorm.cor = 0, mvnorm.mean = 0, mvnorm.cutoff = 0, higher.order.a = 0, higher.order.b = 0),
                         seed = NULL){

  #----------------------------
  # Warning end error messages
  #----------------------------

  if(is.null(N.c)){stop("N.c argument missing, with no default.")}
  if(is.null(Q) & any(sapply(gen.Q, is.null))){stop("Q argument or all gen.Q arguments required.")}
  if(!is.null(Q) & all(!sapply(gen.Q, is.null))){warning("If Q argument is specified, gen.Q arguments are ignored.")}
  if(is.null(IQ)){stop("IQ argument missing, with no default.")}
  if(is.null(VAR)){stop("VAR argument missing, with no default.")}
  if(is.null(gen.att$att.dist)){stop("gen.att$att.dist argument required.")}
  if(is.null(gen.att$mvnorm.cor)){gen.att$mvnorm.cor <- 0}
  if(is.null(gen.att$mvnorm.mean)){gen.att$mvnorm.mean <- 0}
  if(is.null(gen.att$mvnorm.cutoff)){gen.att$mvnorm.cutoff <- 0}
  if(is.null(gen.att$higher.order.a)){gen.att$higher.order.a <- 0}
  if(is.null(gen.att$higher.order.b)){gen.att$higher.order.b <- 0}

  #----------------
  # Gather objects
  #----------------

  if(!is.null(seed)){set.seed(seed)}
  if(is.null(seed)){seed <- sample(1:1000000, size = 1)}

  if(!is.null(Q)){
    arg.Q <- Q
    J <- nrow(Q)
    K <- ncol(Q)
    arg.model <- model
    if(length(model) == 1){model <- rep(model, J)}
  } else {
    arg.Q <- NULL
    J <- gen.Q$J
    K <- gen.Q$K
    Q <- genQ(J, K, gen.Q$nI, gen.Q$propK.J, T, gen.Q$minJ.K, gen.Q$max.Kcor, seed)$Q
    arg.model <- model
    if(length(model) == 1){model <- rep(model, J)}
  }
  if(length(gen.att$mvnorm.mean) == 1){gen.att$mvnorm.mean <- rep(gen.att$mvnorm.mean, K)}
  if(length(gen.att$mvnorm.cutoff) == 1){gen.att$mvnorm.cutoff <- rep(gen.att$mvnorm.cutoff, K)}
  if(length(gen.att$higher.order.a) == 1){gen.att$higher.order.a <- rep(gen.att$higher.order.a, K)}
  if(length(gen.att$higher.order.b) == 1){gen.att$higher.order.b <- rep(gen.att$higher.order.b, K)}
  P <- cbind(P0 = runif(J, (1 - IQ) / 2 - VAR, (1 - IQ) / 2 + VAR),
             P1 = runif(J, (1 - ((1 - IQ) / 2)) - VAR, (1 - ((1 - IQ) / 2)) + VAR))

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

  #---------------
  # Generate data
  #---------------

  sim.v <- NULL
  if(gen.att$att.dist == "uniform"){
    sim.c <- GDINA::simGDINA(N.c, Q, catprob.parm = catprob.parm, att.dist = "uniform")
    if(!is.null(N.v)){
      sim.v <- GDINA::simGDINA(N.v, Q, catprob.parm = catprob.parm, att.dist = "uniform")
    }
  } else if(gen.att$att.dist == "mvnorm"){
    cutoffs <- gen.att$mvnorm.cutoff
    m <- gen.att$mvnorm.mean
    vcov <- matrix(gen.att$mvnorm.cor, K, K); diag(vcov) <- 1
    sim.c <- GDINA::simGDINA(N.c, Q, catprob.parm = catprob.parm, att.dist = "mvnorm",
                             mvnorm.parm = list(mean = m, sigma = vcov, cutoffs = cutoffs))
    if(!is.null(N.v)){
      sim.v <- GDINA::simGDINA(N.v, Q, catprob.parm = catprob.parm, att.dist = "mvnorm",
                               mvnorm.parm = list(mean = m, sigma = vcov, cutoffs = cutoffs))
    }
  } else if(gen.att$att.dist == "higher.order"){
    theta <- rnorm(N.c)
    lambda <- data.frame(a = gen.att$higher.order.a, b = gen.att$higher.order.b)
    sim.c <- GDINA::simGDINA(N.c, Q, catprob.parm = catprob.parm, att.dist = "higher.order",
                             higher.order.parm = list(theta = theta,lambda = lambda))
    if(!is.null(N.v)){
      sim.v <- GDINA::simGDINA(N.v, Q, catprob.parm = catprob.parm, att.dist = "higher.order",
                               higher.order.parm = list(theta = theta,lambda = lambda))
    }
  }

  #----------------
  # Return results
  #----------------

  res <- list()
  res$simdat.c <- sim.c$dat
  res$simalpha.c <- sim.c$attribute
  res$simdat.v <- sim.v$dat
  res$simalpha.v <- sim.v$attribute
  res$simQ <- Q
  res$simcatprob.parm <- catprob.parm
  res$simdelta.parm <- delta.parm
  model[which(rowSums(Q)==1)] <- "GDINA"
  res$specifications <- list(N.c = N.c, N.v = N.v, Q = arg.Q, gen.Q = gen.Q,
                             IQ = IQ, VAR = VAR, model = model, min.param = min.param,
                             gen.att = gen.att, seed = seed)
  return(res)
}
