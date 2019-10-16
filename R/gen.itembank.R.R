#' Item bank generation
#'
#' This function can be used to generate an item bank. The user can provide a Q-matrix or create one defining the number of times each attribute should
#' be measured and the q-vector complexity (e.g., number of attributes that a q-vector can measure). Item parameters are sampled from
#' a uniform distribution with mean = \emph{IQ} and variance = \emph{VAR}.
#'
#' @param Q Q-matrix.
#' @param minJ.K Vector indicating the minimum number of items measuring each attribute.
#' @param complexity Vector indicating the maximum number of attributes being measured by an item in each column of Q. At this moment maximum is 4.
#' @param IQ Item discrimination (mean for the uniform distribution). IQ = \emph{P}(\strong{1}) - \emph{P}(\strong{0}) (Sorrel, Abad, Olea, de la Torre, and Barrada, 2017).
#' @param VAR Item discrimination (variance for the uniform distribution).
#' @param model Vector indicating the model-item correspondnce (0 = one-attribute item, 1 = DINA, 2 = DINO, 3 = A-CDM.
#'
#' @return \code{gen.itembank} returns an object of class \code{gen.Item.Bank}.
#'
#' @references
#'
#' Kaplan, M., de la Torre, J., & Barrada, J. R. (2015). New item selection methods for cognitive diagnosis computerized adaptive testing. \emph{Applied Psychological Measurement, 39}, 167-188.
#'
#' Sorrel, M. A., Abad, F. J., Olea, J., de la Torre, J., & Barrada, J. R. (2017). Inferential item-fit evaluation in cognitive diagnosis modeling. \emph{Applied Psychological Measurement, 41}, 614-631.
#'
#' @examples
#'####################################
#'#        Example 1.                #
#'#     Q and model are provided     #
#'####################################
#' Q <- sim30GDINA$simQ
#' model <- rep(1, each = nrow(Q))
#' IQ <- .70   # P(1), IQ = Low item quality in Kaplan, de la Torre & Barrada (2015)
#' VAR <- 0.10 # High variance in Kaplan et al. (2015)
#' bank <- gen.itembank(Q = Q, IQ = IQ, VAR = VAR, model = model)
#'
#'####################################
#'#        Example 2.                #
#'#   Q and model not are provided   #
#'####################################
#'
#' minJ.K  <- c(50, 50, 50)
#' complexity <- c(3, 3, 3)
#' IQ <- .70   # P(1), IQ = Low item quality in Kaplan, de la Torre & Barrada (2015)
#' VAR <- 0.10 # High variance in Kaplan et al. (2015)
#' bank <- gen.itembank(minJ.K = minJ.K, complexity = complexity, IQ = IQ, VAR = VAR)
#'
#' @export
#'
gen.itembank <- function (Q = NULL, minJ.K = NULL, complexity = NULL, IQ, VAR, model = NULL, ...)
  {

  if((is.null(Q)) & ((is.null(minJ.K)) | (is.null(minJ.K)))) {
    stop("Either Q or minJ.K and complexity are required")
    }

  K1.item.param <- function(IQ, VAR) {
    item.param <- c()
      p1 <- runif(1, IQ - VAR, IQ + VAR)
      p0 <- 1 - p1
      item.param <- c(p0, p1)
    return(item.param)
  }
  DINA.item.param <- function(Q.jj, IQ, VAR) {
    item.param <- c()
      p1 <- runif(1, IQ - VAR, IQ + VAR)
      p0 <- 1 - p1
      if (sum(Q.jj) == 1) {
        item.param <- c(p0, p1)
      }
      if (sum(Q.jj) == 2) {
        item.param <- c(rep(p0, 3), p1)
      }
      if (sum(Q.jj) == 3) {
        item.param <- c(rep(p0, 7), p1)
      }
      if (sum(Q.jj) == 4) {
        item.param <- c(rep(p0, 15), p1)
      }
    return(item.param)
  }
  DINO.item.param <- function(Q.jj, IQ, VAR) {
    item.param <- c()
      p1 <- runif(1, IQ - VAR, IQ + VAR)
      p0 <- 1 - p1
      if (sum(Q.jj) == 1) {
        item.param <- c(p0, p1)
      }
      if (sum(Q.jj) == 2) {
        item.param <- c(p0, rep(p1, 3))
      }
      if (sum(Q.jj) == 3) {
        item.param <- c(p0, rep(p1, 7))
      }
      if (sum(Q.jj) == 4) {
        item.param <- c(p0, rep(p1, 15))
      }
    return(item.param)
  }
  ACDM.item.param <- function(Q.jj, IQ, VAR) {
    item.param <- c()
      p1 <- runif(1, IQ - VAR, IQ + VAR)
      p0 <- 1 - p1
      if (sum(Q.jj) == 1) {
        item.param <- c(p0, p1)
      }
      if (sum(Q.jj) == 2) {
        item.param <- c(p0, rep(p0 + (p1 - p0)/2, 2),
                        p1)
      }
      if (sum(Q.jj) == 3) {
        item.param <- c(p0, rep(p0 + (p1 - p0)/3, 3),
                        rep(p0 + ((p1 - p0)/3) * 2, 3), p1)
      }
      if (sum(Q.jj) == 4) {
        item.param <- c(p0, rep(p0 + (p1 - p0)/4, 4),
                        rep(p0 + ((p1 - p0)/4) * 2, 6), rep(p0 + ((p1 -
                                                                     p0)/4) * 3, 4), p1)
      }
    return(item.param)
  }

  if(is.null(Q)){
  Q <- gen.Q(minJ.K = minJ.K, complexity = complexity)
  J <- nrow(Q)
  J.Kj1 <- which(rowSums(Q) == 1)
  J.Kjc <- which(rowSums(Q) != 1)
  models0 <- sample(x = c(1:3), size = (nrow(Q) - length(J.Kj1)), replace = TRUE) # DINA, DINO, ACDM
  model <- NULL
  model[J.Kj1] <- 0
  model[J.Kjc] <- models0 } else {
    J <- nrow(Q)
    J.Kj1 <- which(rowSums(Q) == 1)
    model[J.Kj1] <- 0
  }

  item.param <- list()
  for (jj in 1:J) {
    Q.jj <- Q[jj, ]
    if (model[[jj]] == 0) {
      item.param[[jj]] <- K1.item.param(IQ, VAR)
    }
    if (model[[jj]] == 1) {
      item.param[[jj]] <- DINA.item.param(Q.jj, IQ, VAR)
    }
    if (model[[jj]] == 2) {
      item.param[[jj]] <- DINO.item.param(Q.jj, IQ, VAR)
    }
    if (model[[jj]] == 3) {
      item.param[[jj]] <- ACDM.item.param(Q.jj, IQ, VAR)
    }
  }

  res <- list("Q" = Q, "item.param" = item.param,
              "model" = model, "options" = list("IQ" = IQ, "VAR" = VAR))

  class(res) <- "gen.Item.Bank"
  return(res)
}
