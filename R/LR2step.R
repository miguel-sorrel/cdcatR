#' Item-level model comparison using 2LR test
#'
#' This function evaluates whether the saturated G-DINA model can be replaced by reduced CDMs without significant loss in model data fit for each item using two-step likelihood ratio test (2LR). Sorrel, de la Torre, Abad, and Olea (2017) and Ma & de la Torre (2018) can be consulted for details.
#'
#' @param fit Calibrated item bank with GDINA or CDM package
#' @param p.adjust.method Correction method for p-values. Possible values include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". See p.adjust function from stats for additional details. Default is holm
#' @param alpha.level Alpha level for decision. Default is 0.05
#' @return \code{LR2step} returns an object of class \code{LR2step}
#' \describe{
#' \item{LR2}{LR2 statistics}
#' \item{pvalues}{p-values associated with the wald statistics}
#' \item{adj.pvalues}{adjusted p-values associated with the wald statistics}
#' \item{df}{degrees of freedom}
#' \item{models.adj.pvalues}{Models selected using the rule \emph{largestp} (Ma, Iaconangelo, & de la Torre, 2016). All statistics whose p-values are less than \code{alpha.level} are rejected. All statistics with p-value larger than \code{alpha.level} define the set of candidate reduced models. The G-DINA model is retained if all statistics are rejected. Whenever the set includes more than one model, the model with the largest p-value was selected as the best model for that item.}
#' }

#' @references
#'
#' Ma, W. & de la Torre, J. (2018). Category-level model selection for the sequential G-DINA model. \emph{Journal of Educational and Behavorial Statistic, 44}, 45-77.
#'
#' Ma, W., Iaconangelo, C., & de la Torre, J. (2016). Model similarity, model selection and attribute classification.
#' \emph{Applied Psychological Measurement, 40}, 200-217.
#'
#' Sorrel, M. A., de la Torre, J., Abad, F. J., & Olea, J. (2017). Two-step likelihood ratio test for item-level model comparison in cognitive diagnosis models. \emph{Methodology, 13}, 39-47.
#'
#' @examples
#' Q <- sim180DINA$simQ
#' dat <- sim180DINA$simdat
#' resGDINA <- GDINA::GDINA(dat = dat, Q = Q, model = "GDINA",verbose = FALSE)
#' resCDM <- CDM::gdina(data = dat, q.matrix = Q, rule = "GDINA", progress = FALSE)
#' LR2.GDINA <- LR_2step(fit = resGDINA)
#' LR2.CDM <- LR_2step(fit = resCDM)
#' table(LR2.GDINA$models.adj.pvalues[which(rowSums(Q) != 1)])
#' table(LR2.CDM$models.adj.pvalues[which(rowSums(Q) != 1)])
#
#' @export
#'

LR_2step <- function(fit, p.adjust.method = "holm", alpha.level = 0.05){
  ## GDINA likelihood is computed from the fit estimates
  # Notes from this function: (ui %*% theta - ci) >= 0

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
