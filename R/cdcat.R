#' @title Cognitively based computerized adaptive test application
#'
#' @description \code{cdcat} conducts a CD-CAT application for a given dataset. Items are selected according to the general discrimination
#' index (GDI; de la Torre & Chiu, 2016; Kaplan, de la Torre, & Barrada, 2016). The next item to be selected by the adaptive algorithm is
#' the one with the highest GDI.
#'
#' @param GDINA.obj Calibrated item bank.
#' @param dat Dataset to be analyzed (if is.null(dat) then dat <- GDINA.obj$options$dat) (i.e., post-hoc CD-CAT simulation).
#' @param MAXJ Maximum number of items to be applied. Default is 20.
#' @param FIXED.LENGTH Fixed CAT-legnth (TRUE) or fixed-precision (FALSE). Default is TRUE.
#' @param att.prior Prior distribution for MAP/EAP estimates.
#' @param post.initial Prior distribution for GDI.
#' @param max.cut Cutoff for fixed-precision (posterior pattern > max.cut). Default is .80.
#' @param i.print Print examinee information. Default is 250.
#'
#' @return \code{cdcat} returns an object of class \code{cdcat}.
#'
#' @references
#'
#' de la Torre, J., & Chiu, C. Y. (2016). General Method of Empirical Q-matrix Validation. \emph{Psychometrika, 81}, 253-273.
#'
#' Kaplan, M., de la Torre, J., & Barrada, J. R. (2015). New item selection methods for cognitive diagnosis computerized adaptive testing. \emph{Applied Psychological Measurement, 39}, 167-188.
#'
#' @examples
#'####################################
#'#        Example 1.                #
#'#     CD-CAT simulation for a      #
#'#     GDINA obj                    #
#'####################################
#'
#'#--------Data generation --------#
#'Q <- sim180GDINA$simQ
#'K <- ncol(Q)
#'dat <- sim180GDINA$simdat
#'att <- sim180GDINA$simalpha
#'
#'#--------Model estimation --------#
#'fit <- GDINA(dat = dat, Q = Q, verbose = 0)
#'
#'#--------CD-CAT --------#
#'res.FIXJ <- cdcat(dat = fit$options$dat,
#'                  GDINA.obj = fit,
#'                  FIXED.LENGTH = TRUE)
#'res.VARJ <- cdcat(dat = fit$options$dat,
#'                  GDINA.obj = fit,
#'                  FIXED.LENGTH = FALSE)
#'
#'#--------Results --------#
#'res.FIXJ$est[[1]] # estimates for the first examinee (fixed-length)
#'res.VARJ$est[[1]] # estimates for the first examinee (fixed-precision)
#'att.plot(res.FIXJ, i = 1) # plot for estimates for the first examinee (fixed-length)
#'att.plot(res.VARJ, i = 1) # plot for estimates for the first examinee (fixed-length)
#'# FIXJ summary
#'res.FIXJ.sum.real <- cdcat.summary(cdcat.obj = res.FIXJ,
#'                                   alpha = att) # vs. real accuracy
#'res.FIXJ.sum.real$recovery$plotPCV
#'res.FIXJ.sum.real$recovery$plotPCA
#'# VARJ summary
#'res.VARJ.sum.post <- cdcat.summary(cdcat.obj = res.VARJ, alpha = att)
#'res.VARJ.sum.post$CATlength$stats
#'res.VARJ.sum.post$CATlength$plot
#'res.VARJ.sum.post$recovery
#'# Post-hoc CAT simulation (only if dat is fit$options$dat)
#'att.J <- personparm(fit, "MAP")[, -(K+1)]
#'class.J <- ClassRate(att, att.J) # upper-limit for accuracy
#'res.FIXJ.sum.post <- cdcat.summary(cdcat.obj = res.FIXJ, alpha = att.J)
#'res.FIXJ.sum.post$recovery$plotPCV + geom_hline(yintercept = class.J$PCV[K], color = "red")
#'res.FIXJ.sum.post$recovery$plotPCA + geom_hline(yintercept = class.J$PCA, color = "red")
#'
#'####################################
#'# Example 2. #
#'# CD-CAT simulation for #
#'# multiple GDINA objs and #
#'# comparison of performance on a
#'# validation sample #
#'####################################
#'
#'#--------Data --------#
#'Q <- sim155complex$simQ
#'K <- ncol(Q)
#'parm <- sim155complex$simcatprob.parm
#'dat.c <- sim155complex$simdat.c
#'att.c <- sim155complex$simalpha.c
#'dat.v <- sim155complex$simdat.v
#'att.v <- sim155complex$simalpha.v
#'
#'#--------(multiple) Model estimation --------#
#'fitTRUE <- GDINA(dat = dat.c, Q = Q, catprob.parm = parm, control = list(maxitr = 0), verbose = 0)
#'fitGDINA <- GDINA(dat = dat.c, Q = Q, verbose = 0)
#'fitDINA <- GDINA(dat = dat.c, Q = Q, model = "DINA", verbose = 0)
#'fitDINO <- GDINA(dat = dat.c, Q = Q, model = "DINO", verbose = 0)
#'fitACDM <- GDINA(dat = dat.c, Q = Q, model = "ACDM", verbose = 0)
#'LR2 <- modelcomp(GDINA.obj = fitGDINA, method = "LR",
#'                 LR.args = list(LR.approx = TRUE), models = c("DINA", "DINO","ACDM"))
#'alpha.level <- 0.05
#'model <- apply(LR2$pvalues, 1, function(x) {
#'  if (max(x, na.rm = TRUE) > alpha.level/(sum(rowSums(Q)>1)*3)) {
#'    which.max(x)}
#'  else {
#'    return(0)}})
#'models <- rep(0, nrow(Q))
#'models[which(rowSums(Q) != 1)] <- model
#'models <- models
#'fitLR2 <- GDINA(dat = dat.c, Q = Q, model = models, verbose = 0)
#'
#'#--------CD-CAT --------#
#'fit.l <- list(fitTRUE, fitGDINA, fitDINA, fitDINO, fitACDM, fitLR2)
#'res.FIXJ.l <- res.VARJ.l <- list()
#'for(mm in 1:length(fit.l)) {
#'  fit <- fit.l[[mm]]
#'  res.FIXJ.l[[mm]] <- cdcat(dat = dat.v,
#'                            GDINA.obj = fit,
#'                          FIXED.LENGTH = TRUE)
#'  res.VARJ.l[[mm]] <- cdcat(dat = dat.v,
#'                          GDINA.obj = fit,
#'                            FIXED.LENGTH = FALSE)
#'}
#'
#'#--------Results --------#
#'fitbest <- GDINA(dat = dat.v, Q = Q, catprob.parm = parm, control = list(maxitr = 1), verbose = 0)
#'fitbest.acc <- personparm(fitbest, "MAP")[, -(K+1)]
#'class.J <- ClassRate(att.v, fitbest.acc) # upper-limit for accuracy
#'# FIXJ comparison
#'res.FIXJ.sum.post.comp <- cdcat.comp(cdcat.obj.l = res.FIXJ.l, alpha = att.v)
#'res.FIXJ.sum.post.comp$PCVcomp + geom_hline(yintercept = class.J$PCV[K], color = "red")
#'res.FIXJ.sum.post.comp$PCAmcomp + geom_hline(yintercept = class.J$PCA, color = "red")
#'# VARJ comparison
#'res.VARJ.sum.post.comp <- cdcat.comp(cdcat.obj.l = res.VARJ.l, alpha = att.v)
#'res.VARJ.sum.post.comp$stats
#'res.VARJ.sum.post.comp$plots
#'res.VARJ.sum.post.comp$recovery
#'
#' @export
#'
cdcat <- function(GDINA.obj,
                  dat = NULL,
                  MAXJ = 20,
                  FIXED.LENGTH = TRUE,
                  att.prior = NULL,
                  post.initial = NULL,
                  max.cut = 0.80,
                  i.print = 250, ...)
  {

  GDI.M <- function(LC.prob, prior) {
    GDI.est <- c(0)
    for (jj in 1:nrow(LC.prob)) {
      MPj <- sum(LC.prob[jj,]*prior)
      GDI.est[jj] <- sum(prior*(LC.prob[jj,]-MPj)^2)    }
    return(GDI.est)
  }

  est <- GDINA.obj
  if(length(table(est$model)) == 1) {model <- names(table(est$model))} else {model <- "Combination"}

  if (is.null(dat)) {
    X <- GDINA.obj$options$dat
  } else {
    X <- dat
  }

  Q <- GDINA.obj$options$Q
  N <- nrow(X)
  K <- ncol(Q)
  L <- 2^K
  pattern <- attributepattern(K)
  Lclass <- apply(attributepattern(K),MARGIN = 1,FUN =  function(x) {paste(x,collapse = "")})

  if (is.null(att.prior)) {att.prior <- rep(1/L, L)}
  if (is.null(post.initial)) {post.initial <- matrix(data = rep(rep(1/L, L), N), nrow = N)}

  out <- list()

  for(i in 1:N) {
    if (i %% i.print == 0) {print(i)}

    mlogPost_GDI <- post.initial[i,]
    est.GDI <- NULL
    item.log <- NULL
    GDI <- GDI.M(LC.prob = est$LC.prob,  mlogPost_GDI)
    names(GDI) <- 1:length((GDI))
    jjselect <- 1
    jjcatGDI <- IP <- NULL

    if (FIXED.LENGTH == TRUE) {
      for (jjselect in 1:MAXJ) {

        if (is.null(jjcatGDI)) {
          jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI)))
        } else {
          jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI[-c(jjcatGDI)])))
        }

        x.jj = X[i,jjcatGDI[jjselect]]
        IP = cbind(IP,x.jj)

        p.xi <- matrix(est$LC.prob[jjcatGDI,], ncol = L)
        p.xi.aj1_GDI <- matrix(0,nrow = L, ncol = length(IP))
        for (l in 1L:L){
          for (jj in 1L:length(IP)){
            res <- 1
            p.xi.aj1_GDI[l,jj] = res*(p.xi[jj,l]^IP[jj])*((1-p.xi[jj,l])^(1-IP[jj]))
            res <- p.xi.aj1_GDI[l,jj]
          }
        }

        p.xi.aj_GDI <- (t(as.matrix(apply(p.xi.aj1_GDI, 1, prod))))
        mPost_GDI   <- p.xi.aj_GDI * att.prior
        msumPost_GDI <- sum(mPost_GDI)
        mlogPost_GDI <- mPost_GDI/msumPost_GDI
        p.aj_GDI <- mlogPost_GDI %*% pattern

        GDI <- GDI.M(LC.prob = est$LC.prob, prior = mlogPost_GDI)
        names(GDI) <- 1:length((GDI))

        est.GDI <- rbind(est.GDI,
                         c(Lclass[which.max(p.xi.aj_GDI)], max(p.xi.aj_GDI),
                           Lclass[which.max(mlogPost_GDI)], max(mlogPost_GDI),
                           paste(as.numeric(p.aj_GDI > .50), collapse = ""), p.aj_GDI))

        jjselect <- jjselect + 1
        if((jjselect) > MAXJ) {break}
      }

      colnames(est.GDI) <- c("ML", "Lik", "MAP", "Post", "EAP", paste("K", 1:K, sep = ""))

    } else {

      while(max(mlogPost_GDI) < max.cut){

        if (is.null(jjcatGDI)) {
          jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI)))
        } else {
          jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI[-c(jjcatGDI)])))
        }

        x.jj = X[i,jjcatGDI[jjselect]]
        IP = cbind(IP,x.jj)

        p.xi <- matrix(est$LC.prob[jjcatGDI,], ncol = L)
        p.xi.aj1_GDI <- matrix(0,nrow = L, ncol = length(IP))
        for (l in 1L:L){
          for (jj in 1L:length(IP)){
            res <- 1
            p.xi.aj1_GDI[l,jj] = res*(p.xi[jj,l]^IP[jj])*((1-p.xi[jj,l])^(1-IP[jj]))
            res <- p.xi.aj1_GDI[l,jj]
          }
        }

        p.xi.aj_GDI <- (t(as.matrix(apply(p.xi.aj1_GDI, 1, prod))))
        mPost_GDI   <- p.xi.aj_GDI * att.prior
        msumPost_GDI <- sum(mPost_GDI)
        mlogPost_GDI <- mPost_GDI/msumPost_GDI
        p.aj_GDI <- mlogPost_GDI %*% pattern

        GDI <- GDI.M(LC.prob = est$LC.prob, prior = mlogPost_GDI)
        names(GDI) <- 1:length((GDI))

        est.GDI <- rbind(est.GDI,
                         c(Lclass[which.max(p.xi.aj_GDI)], max(p.xi.aj_GDI),
                           Lclass[which.max(mlogPost_GDI)], max(mlogPost_GDI),
                           paste(as.numeric(p.aj_GDI > .50), collapse = ""), p.aj_GDI))
        jjselect <- jjselect + 1
        if((jjselect) > MAXJ) {break}
      }

      colnames(est.GDI) <- c("ML", "Lik", "MAP", "Post", "EAP", paste("K", 1:K, sep = ""))
    }


    item.log <- c(item.log, jjcatGDI)
    out[[i]] <- list(est.cat.GDI = est.GDI,
                     item.usage = item.log)

  }
  res <- list()
  res$est <- out
  res$specifications <- list("dat" = dat, "GDINA.obj" = GDINA.obj, "model" = model, "MAXJ" = MAXJ, "FIXED.LENGTH" = FIXED.LENGTH,
                             "att.prior" = att.prior, "post.initial" = post.initial, "max.cut" = max.cut,
                             "i.print" = i.print)
  class(res) <- "cdcat"

  return(res)
}
