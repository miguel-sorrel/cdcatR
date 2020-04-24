#' Item-level model comparison using 2LR test
#'
#' This function evaluates whether the saturated G-DINA model can be replaced by reduced CDMs without significant loss in model data fit for each item using two-step likelihood ratio test (2LR). Sorrel, de la Torre, Abad, and Olea (2017) and Ma & de la Torre (2018) can be consulted for details.
#'
#' @param fit Calibrated item bank with the \code{GDINA::GDINA} (Ma & de la Torre, 2020) or \code{CDM::gdina} (Robitzsch et al., 2020) R packages functions
#' @param p.adjust.method Character vector of length 1. Correction method for \emph{p}-values. Possible values include \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}, and \code{"none"}. See \code{p.adjust} function from the \code{stats} R package for additional details. Default is \code{holm}
#' @param alpha.level Numeric vector of length 1. Alpha level for decision. Default is \code{0.05}
#'
#' @return \code{LR2.step} returns an object of class \code{LR2.step}
#' \describe{
#' \item{LR2}{Numeric matrix. LR2 statistics}
#' \item{pvalues}{Numeric matrix. \emph{p}-values associated with the 2LR statistics}
#' \item{adj.pvalues}{Numeric matrix. Adjusted \emph{p}-values associated with the 2LR statistics}
#' \item{df}{Numeric matrix. Degrees of freedom}
#' \item{models.adj.pvalues}{Character vector denoting the model selected for each item using the \emph{largestp} rule (Ma et al., 2016). All statistics whose \emph{p}-values are less than \code{alpha.level} are rejected. All statistics with \emph{p}-value larger than \code{alpha.level} define the set of candidate reduced models. The G-DINA model is retained if all statistics are rejected. Whenever the set includes more than one model, the model with the largest \emph{p}-value was selected as the best model for that item}
#' }

#' @references
#'
#' Ma, W. & de la Torre, J. (2018). Category-level model selection for the sequential G-DINA model. \emph{Journal of Educational and Behavorial Statistic, 44}, 45-77.
#'
#' Ma, W. & de la Torre, J. (2020). GDINA: The generalized DINA model framework. R package version 2.7.9. Retrived from https://CRAN.R-project.org/package=GDINA
#'
#' Ma, W., Iaconangelo, C., & de la Torre, J. (2016). Model similarity, model selection and attribute classification.
#' \emph{Applied Psychological Measurement, 40}, 200-217.
#'
#' Robitzsch, A., Kiefer, T., George, A. C., & Uenlue, A. (2020). CDM: Cognitive Diagnosis Modeling. R package version 7.5-15. https://CRAN.R-project.org/package=CDM
#'
#' Sorrel, M. A., de la Torre, J., Abad, F. J., & Olea, J. (2017). Two-step likelihood ratio test for item-level model comparison in cognitive diagnosis models. \emph{Methodology, 13}, 39-47.
#'
#' @import stats
#'
#' @examples
#'Q <- sim180DINA$simQ
#'dat <- sim180DINA$simdat
#'resGDINA <- GDINA::GDINA(dat = dat, Q = Q, model = "GDINA",verbose = FALSE)
#'resCDM <- CDM::gdina(data = dat, q.matrix = Q, rule = "GDINA", progress = FALSE)
#'LR2.GDINA <- LR.2step(fit = resGDINA)
#'LR2.CDM <- LR.2step(fit = resCDM)
#'mean(LR2.GDINA$models.adj.pvalues[which(rowSums(Q) != 1)] ==
#'       sim180DINA$specifications$item.bank$specifications$model[which(rowSums(Q) != 1)])
#'mean(LR2.CDM$models.adj.pvalues[which(rowSums(Q) != 1)] ==
#'       sim180DINA$specifications$item.bank$specifications$model[which(rowSums(Q) != 1)])
#'
#' @export
#'

LR.2step <- function(fit, p.adjust.method = "holm", alpha.level = 0.05)
  {

  make.Lik.DINA.k <- function(vP){
    p <- rep(vP[1], classes - 1)
    p[classes] <- vP[1] + vP[2]
    logp <- log(p)
    logq <- log(1 - p)
    L <- sum(RNjjob[1,] * logp) + sum(RNjjob[3,] * logq)
    return(-L)
  }
  make.Lik.DINO.k <- function(vP){
    p <- rep(vP[1] + vP[2], classes)
    p[1] <- vP[1]
    logp <- log(p)
    logq <- log(1 - p)
    L <- sum(RNjjob[1,] * logp) + sum(RNjjob[3,] * logq)
    return(-L)
  }
  make.Lik.ACDM.k <- function(vP){
    par <- matrix(vP[-1], ncol = 1)
    p <- c(vP[1] + posible.patterns %*% par)
    logp <- log(p)
    logq <- log(1 - p)
    L <- sum(RNjjob[1,] * logp) + sum(RNjjob[3,] * logq)
    return(-L)
  }
  make.Lik.GDINA.k <- function(pars){
    p <- pars
    logp <- log(p)
    logq <- log(1 - p)
    L <- sum(RNjjob[1,] * logp) + sum(RNjjob[3,] * logq)
    return(-L)
  }

  est <- fit
  if(!is.null(est)){
    if(!is.null(est$extra$call)){ # package GDINA
      Q <- est$options$Q
      attpat <- GDINA::attributepattern(ncol(Q))
      Rlj <- est$technicals$expectedCorrect
      Nlj <- est$technicals$expectedTotal
      catprob.parm <- est$catprob.parm
    } else { # package CDM
      Q <- est$q.matrix
      attpat <- GDINA::attributepattern(ncol(Q))
      attpat.s <- apply(attpat, 1, paste, collapse = "")
      RLj <- est$control$R.lj
      RLj <- RLj[,match(attpat.s, colnames(RLj))]
      NLj <- est$control$I.lj
      NLj <- NLj[,match(attpat.s, colnames(RLj))]
      Rlj <- Nlj <- matrix(NA, nrow = nrow(Q), ncol = 2^max(rowSums(Q)))
      catprob.parm <- list()
      for(j in 1:nrow(Q)){
        kj <- which(Q[j,] == 1)
        tmp <- factor(apply(attpat[,kj, drop = F], 1, paste, collapse = ""))
        tmp <- factor(tmp, levels = apply(unique(GDINA::attributepattern(length(kj))), 1, paste, collapse = ""))
        for(l in 1:(2^length(kj))){
          Rlj[j,l] <- sum(RLj[j, which(as.numeric(tmp) == l)])
          Nlj[j,l] <- sum(NLj[j, which(as.numeric(tmp) == l)])
        }
        catprob.parm[[j]] <- est$probitem[est$probitem$itemno == j,]$prob
      }
      Rlj[is.na(Rlj)] <- 0
      Nlj[is.na(Nlj)] <- 0
    }
  }

  Kjs <- rowSums(Q)
  Ks <- cumsum(2^Kjs)
  item <- which(rowSums(Q) > 1)

  J <- nrow(Q)
  out <- list()
  LR2 <- matrix(data = 0, nrow = J, ncol = 4)
  LR2.p <- matrix(data = 0, nrow = J, ncol = 4)
  df <- matrix(data = 0, nrow = J, ncol = 4)

  for(jj in 1:J){
    pars.jj <- catprob.parm[[jj]]
    pars.jj[pars.jj < 0] <- 0
    pars.jj[pars.jj > 1] <- 1
    Kjj <- sum(Q[jj, ])
    Kj <- length(pars.jj)
    RNjjob <- rbind(Rlj[jj, 1:Kj], Nlj[jj, 1:Kj], Nlj[jj, 1:Kj] - Rlj[jj, 1:Kj])
    classes <- 2^Kjs[jj]
    posible.patterns <- GDINA::attributepattern(Kjj)

    par.DINA <- constrOptim(theta = c(0.1, 0.1), f = make.Lik.DINA.k, grad = NULL,
                            ui = rbind(c(1, 0), c(0, 1), c(1, 1), c(-1, 0), c(0, -1), c(-1, -1)),
                            ci = c(.0001, .0001, .0001, -.9999, -.9999, -.9999))

    par.DINO <- constrOptim(theta = c(0.1, 0.1), f = make.Lik.DINO.k, grad = NULL,
                            ui = rbind(c(1, 0), c(0, 1), c(1, 1), c(-1, 0), c(0, -1), c(-1, -1)),
                            ci = c(.0001, .0001, .0001, -.9999, -.9999 , -.9999))

    par.ACDM <- constrOptim(theta = rep(0.1, Kjs[jj]+1), f = make.Lik.ACDM.k, grad = NULL,
                            ui = rbind(diag(Kjs[jj] + 1), -diag(Kjs[jj] + 1), rep(-1, Kjs[jj] + 1)),
                            ci = c(rep(.0001, Kjs[jj] + 1), rep(-.9999, Kjs[jj] + 2)))

    par.GDINA.value <- make.Lik.GDINA.k(pars = pars.jj)

    if(ncol(RNjjob) == 2){
      LR2[jj,] <- LR2.p[jj,] <- df[jj,] <- cbind(jj, NA, NA, NA)
    } else {
      LR2DINA <- 2 * (par.GDINA.value * (-1) - par.DINA$value * (-1))
      LR2DINO <- 2 * (par.GDINA.value * (-1) - par.DINO$value * (-1))
      LR2ACDM <- 2 * (par.GDINA.value * (-1) - par.ACDM$value * (-1))

      dfDINA  <- dfDINO <- ncol(RNjjob) - 2
      dfACDM  <- ncol(RNjjob) - (Kjj + 1)

      pDINA <-  1 - pchisq(LR2DINA, dfDINA)
      pDINO <-  1 - pchisq(LR2DINO, dfDINO)
      pACDM <-  1 - pchisq(LR2ACDM, dfACDM)

      LR2[jj,] <- cbind(jj, LR2DINA, LR2DINO, LR2ACDM)
      LR2.p[jj,] <- cbind(jj, pDINA, pDINO, pACDM)
      df[jj,] <- cbind(jj, dfDINA, dfDINO, dfACDM)
    }
  }

  if((sum(na.omit(LR2[, -1] == -Inf)) > 0) | (sum(is.nan(LR2[, -1])) > 0)) {
  prob.items <- apply(LR2[, -1], 2, function(x) {which(x == -Inf)})[, 1]
    warning(c("2LR statistic couldn't be computed for items ", paste(prob.items, collapse = ", "),
              ". G-DINA was retained for those items"))
  }

  LR2[, -1][which(LR2[, -1] == -Inf)] <- -9999
  LR2[, -1][is.nan(LR2[, -1])] <- -9999
  LR2.p[, -1][is.nan(LR2.p[, -1])] <- 0

  LR2 <- na.omit(LR2)[, -1]
  LR2.p <- na.omit(LR2.p)[, -1]
  df <- na.omit(df)[, -1]
  LR2.adjp <- LR2.p
  LR2.adjp <- matrix(p.adjust(LR2.adjp, method = p.adjust.method), ncol = ncol(LR2.p), nrow = nrow(LR2.p))

  model.alpha <- apply(LR2.adjp, 1, function(x){if(max(x, na.rm = T) > alpha.level){NA} else {return(0)}})
  model.alpha[is.na(model.alpha)] <- apply(LR2.p[is.na(model.alpha),], 1, which.max)

  models <- rep(0, nrow(Q))
  models[which(rowSums(Q) != 1)] <- model.alpha
  models.adjp <- models

  models.adjp[which(models.adjp == 0)] <- "GDINA"
  models.adjp[which(models.adjp == 1)] <- "DINA"
  models.adjp[which(models.adjp == 2)] <- "DINO"
  models.adjp[which(models.adjp == 3)] <- "ACDM"

  colnames(LR2) <- colnames(LR2.p) <- colnames(LR2.adjp) <- colnames(df) <- c("DINA", "DINO", "ACDM")
  rownames(LR2) <- rownames(LR2.p) <- rownames(LR2.adjp) <- rownames(df) <- paste("Item", item)
  out <- list(LR2 = LR2, pvalues = LR2.p, adj.pvalues = LR2.adjp, df = df, models.adj.pvalues = models.adjp)

  return(out)
}
