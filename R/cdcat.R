#' @title Cognitively based computerized adaptive test application
#'
#' @description \code{cdcat} conducts a CD-CAT application for a given dataset. Different item selection rules can be used: the general discrimination index
#' (GDI; de la Torre & Chiu, 2016; Kaplan et al., 2015), the Jensen-Shannon divergence index (JSD; Kang et al., 2017; Minchen & de la Torre, 2016;
#' Yigit et al., 2018), the posterior-weighted Kullback-Leibler index (PWKL; Cheng, 2009), the modified PWKL index (MPWKL; Kaplan et al., 2015), the nonparametric
#' item selection method (NPS; Chang et al., 2018), or random selection. Fixed length or fixed precision CD-CAT can be applied. Fixed precision CD-CAT with
#' NPS is available, by using the pseudo-posterior probability of each student mastering each attribute (experimental).
#'
#' @param fit Calibrated item bank with the GDINA or CDM package
#' @param dat Dataset to be analyzed. If \code{is.null(dat)} it takes data from the fit object (i.e., post-hoc simulation)
#' @param Q Q-matrix to be used in the analysis. If \code{is.null(Q)} it takes Q from the fit object
#' @param itemSelect Item selection rule: GDI, JSD, MPWKL, PWKL, NPS, random
#' @param MAXJ Maximum number of items to be applied. Default is 20
#' @param FIXED.LENGTH Fixed CAT-length (TRUE) or fixed-precision (FALSE). Default is TRUE
#' @param att.prior Prior distribution for MAP/EAP estimates
#' @param post.initial Prior distribution for \code{itemSelect}
#' @param max.cut Cutoff for fixed-precision (posterior pattern > max.cut). When \code{itemSelect = "NPS"} this is evaluated at the attribute level using the pseudo-posterior probabilities for each attribute (\emph{K} posterior probabilities > max.cut). Default is .80. A higher cutoff is recommended when \code{itemSelect = "NPS"}
#' @param NPS.args A list of options when \code{itemSelect = "NPS"}. \code{gate} = "AND" or "OR", depending on whether a conjunctive o disjunctive nonparametric CDM is used. \code{pseudo.prob} = pseudo-posterior probability of each examinee mastering each attribute (experimental). \code{w.type} = weight type used for computing the pseudo-posterior probability (experimental). \code{seed} = NPS has a random component, so a seed is required for consistent results.
#' @param n.cores Number of cores to be used during parallelization. Default is 2
#' @param i.print Print examinee information. Default is 250
#'
#' @return \code{cdcat} returns an object of class \code{cdcat}.
#' \describe{
#' \item{est}{A list of that contains for each examinee the mastery probability estimates at each step of the CAT \code{est.cat} and the items applied \code{item.usage}}
#' \item{specifications}{A list of that contains all the specifications}
#' }

#' @references
#'
#' Chang, Y.-P., Chiu, C.-Y., & Tsai, R.-C. (2019). Nonparametric CAT for CD in educational settings with small samples. \emph{Applied Psychological Measurement, 43}, 543-561.
#'
#' Cheng, Y. (2009). When cognitive diagnosis meets computerized adaptive testing: CD-CAT. \emph{Psychometrika, 74}, 619-632.
#'
#' de la Torre, J., & Chiu, C. Y. (2016). General method of empirical Q-matrix validation. \emph{Psychometrika, 81}, 253-273.
#'
#' Kang, H.-A., Zhang, S., & Chang, H.-H. (2017). Dual-objective item selection criteria in cognitive diagnostic computerized adaptive testing. \emph{Journal of Educational Measurement, 54}, 165-183.
#'
#' Kaplan, M., de la Torre, J., & Barrada, J. R. (2015). New item selection methods for cognitive diagnosis computerized adaptive testing. \emph{Applied Psychological Measurement, 39}, 167-188.
#'
#' Minchen, N., & de la Torre, J. (2016, July). \emph{The continuous G-DINA model and the Jensen-Shannon divergence}. Paper presented at the International Meeting of the Psychometric Society, Asheville, NC, United States.
#'
#' Yigit, H. D., Sorrel, M. A., de la Torre, J. (2018). Computerized adaptive testing for cognitively based multiple-choice data. \emph{Applied Psychological Measurement, 43}, 388-401.
#'
#' @import foreach
#' @import stats
#'
#' @examples
#'######################################
#'# Example 1.                         #
#'# CD-CAT simulation for a GDINA obj  #
#'######################################
#'
#'#-----------Data generation----------#
#'Q <- sim180GDINA$simQ
#'K <- ncol(Q)
#'dat <- sim180GDINA$simdat
#'att <- sim180GDINA$simalpha
#'
#'#----------Model estimation----------#
#'fit <- GDINA::GDINA(dat = dat, Q = Q, verbose = 0) # GDINA package
#'#fit <- CDM::gdina(data = dat, q.matrix = Q, progress = 0) # CDM package
#'
#'#---------------CD-CAT---------------#
#'res.FIXJ <- cdcat(fit = fit, dat = dat, FIXED.LENGTH = TRUE, n.cores = 4)
#'res.VARJ <- cdcat(fit = fit, dat = dat, FIXED.LENGTH = FALSE, n.cores = 4)
#'
#'#---------------Results--------------#
#'res.FIXJ$est[[1]] # estimates for the first examinee (fixed-length)
#'res.VARJ$est[[1]] # estimates for the first examinee (fixed-precision)
#'att.plot(cdcat.obj = res.FIXJ, i = 1) # plot for the first examinee (fixed-length)
#'att.plot(cdcat.obj = res.VARJ, i = 1) # plot  for the first examinee (fixed-precision)
#'# FIXJ summary
#'res.FIXJ.sum.real <- cdcat.summary(cdcat.obj = res.FIXJ, alpha = att) # vs. real accuracy
#'res.FIXJ.sum.real$recovery$plotPCV
#'res.FIXJ.sum.real$recovery$plotPCA
#'# VARJ summary
#'res.VARJ.sum.post <- cdcat.summary(cdcat.obj = res.VARJ, alpha = att)
#'res.VARJ.sum.post$CATlength$stats
#'res.VARJ.sum.post$CATlength$plot
#'res.VARJ.sum.post$recovery
#'# Post-hoc CAT simulation (only if dat is fit$options$dat)
#'att.J <- GDINA::personparm(fit, "MAP")[, -(K+1)] # GDINA package
#'# att.J <- t(sapply(strsplit(as.character(fit$pattern$map.est), ""), as.numeric)) # CDM package
#'class.J <- GDINA::ClassRate(att, att.J) # upper-limit for accuracy
#'res.FIXJ.sum.post <- cdcat.summary(cdcat.obj = res.FIXJ, alpha = att.J)
#'res.FIXJ.sum.post$recovery$plotPCV + ggplot2::geom_hline(yintercept = class.J$PCV[K], color = "red")
#'res.FIXJ.sum.post$recovery$plotPCA + ggplot2::geom_hline(yintercept = class.J$PCA, color = "red")
#'
#'######################################
#'# Example 2.                         #
#'# CD-CAT simulation for multiple     #
#'# GDINA objs and comparison of       #
#'# performance on a validation sample #
#'######################################
#'
#'#----------------Data----------------#
#'Q <- sim155complex$simQ
#'K <- ncol(Q)
#'parm <- sim155complex$simcatprob.parm
#'dat.c <- sim155complex$simdat.c
#'att.c <- sim155complex$simalpha.c
#'dat.v <- sim155complex$simdat.v
#'att.v <- sim155complex$simalpha.v
#'
#'#-----(multiple) Model estimation----#
#'fitTRUE <- GDINA::GDINA(dat = dat.c, Q = Q, catprob.parm = parm,
#'
#'            control = list(maxitr = 0), verbose = 0)
#'fitGDINA <- GDINA::GDINA(dat = dat.c, Q = Q, verbose = 0)
#'fitDINA <- GDINA::GDINA(dat = dat.c, Q = Q, model = "DINA", verbose = 0)
#'fitDINO <- GDINA::GDINA(dat = dat.c, Q = Q, model = "DINO", verbose = 0)
#'fitACDM <- GDINA::GDINA(dat = dat.c, Q = Q, model = "ACDM", verbose = 0)
#'LR2step <- LR.2step(fitGDINA)
#'models <- LR2step$models.adj.pvalues
#'fitLR2 <- GDINA::GDINA(dat = dat.c, Q = Q, model = models, verbose = 0)
#'
#'#---------------CD-CAT---------------#
#'fit.l <- list(fitTRUE, fitGDINA, fitDINA, fitDINO, fitACDM, fitLR2)
#'res.FIXJ.l <- res.VARJ.l <- list()
#'for(mm in 1:length(fit.l)) {
#'  fit <- fit.l[[mm]]
#'  res.FIXJ.l[[mm]] <- cdcat(fit = fit, dat = dat.v, FIXED.LENGTH = TRUE, n.cores = 4)
#'  res.VARJ.l[[mm]] <- cdcat(fit = fit, dat = dat.v, FIXED.LENGTH = FALSE, n.cores = 4)
#'}
#'
#'#---------------Results--------------#
#'fitbest <- GDINA::GDINA(dat = dat.v, Q = Q, catprob.parm = parm,
#'           control = list(maxitr = 1), verbose = 0)
#'fitbest.acc <- GDINA::personparm(fitbest, "MAP")[, -(K+1)]
#'class.J <- GDINA::ClassRate(att.v, fitbest.acc) # upper-limit for accuracy
#'# FIXJ comparison
#'res.FIXJ.sum.post.comp <- cdcat.comp(cdcat.obj.l = res.FIXJ.l, alpha = att.v)
#'res.FIXJ.sum.post.comp$PCVcomp + ggplot2::geom_hline(yintercept = class.J$PCV[K], color = "red")
#'res.FIXJ.sum.post.comp$PCAmcomp + ggplot2::geom_hline(yintercept = class.J$PCA, color = "red")
#'# VARJ comparison
#'res.VARJ.sum.post.comp <- cdcat.comp(cdcat.obj.l = res.VARJ.l, alpha = att.v)
#'res.VARJ.sum.post.comp$stats
#'res.VARJ.sum.post.comp$plots
#'res.VARJ.sum.post.comp$recovery
#'
#'######################################
#'# Example 3.                         #
#'# Nonparametric CD-CAT for           #
#'# small-scale assessment             #
#'######################################
#'
#'#-----------Data generation----------#
#'Q <- sim180DINA$simQ
#'K <- ncol(Q)
#'N <- 50
#'dat <- sim180DINA$simdat[1:N,]
#'att <- sim180DINA$simalpha[1:N,]
#'
#'#--------Nonparametric CD-CAT--------#
#'res.NPS.FIXJ <- cdcat(dat = dat, Q = Q, itemSelect = "NPS", FIXED.LENGTH = TRUE,
#'                      NPS.args = list(gate = "AND", pseudo.prob = TRUE, w.type = 1, seed = 12345),
#'                      n.cores = 4)
#'res.NPS.VARJ <- cdcat(dat = dat, Q = Q, itemSelect = "NPS", FIXED.LENGTH = FALSE, max.cut = 0.95,
#'                      NPS.args = list(gate = "AND", pseudo.prob = TRUE, w.type = 1, seed = 12345),
#'                      n.cores = 4)
#'
#'#---------------Results--------------#
#'res.NPS.FIXJ$est[[1]] # estimates for the first examinee (fixed-length)
#'res.NPS.VARJ$est[[1]] # estimates for the first examinee (fixed-precision)
#'att.plot(res.NPS.FIXJ, i = 1) # plot for estimates for the first examinee (fixed-length)
#'att.plot(res.NPS.VARJ, i = 1) # plot for estimates for the first examinee (fixed-precision)
#'# FIXJ summary
#'res.NPS.FIXJ.sum.real <- cdcat.summary(cdcat.obj = res.NPS.FIXJ, alpha = att) # vs. real accuracy
#'res.NPS.FIXJ.sum.real$recovery$plotPCV
#'res.NPS.FIXJ.sum.real$recovery$plotPCA
#'# VARJ summary
#'res.NPS.VARJ.sum.post <- cdcat.summary(cdcat.obj = res.NPS.VARJ, alpha = att)
#'res.NPS.VARJ.sum.post$CATlength$stats
#'res.NPS.VARJ.sum.post$CATlength$plot
#'res.NPS.VARJ.sum.post$recovery
#'# Post-hoc CAT simulation
#'fit <- NPCD::AlphaNP(Y = dat, Q = Q, gate = "AND")
#'att.J <- fit$alpha.est
#'class.J <- GDINA::ClassRate(att, att.J) # upper-limit for accuracy
#'res.NPS.FIXJ.sum.post <- cdcat.summary(cdcat.obj = res.NPS.FIXJ, alpha = att.J)
#'res.NPS.FIXJ.sum.post$recovery$plotPCV + ggplot2::geom_hline(yintercept = class.J$PCV[K],
#'                                                             color = "firebrick3")
#'res.NPS.FIXJ.sum.post$recovery$plotPCA + ggplot2::geom_hline(yintercept = class.J$PCA,
#'                                                             color = "firebrick3")
#'
#' @export
#'
cdcat <- function(fit = NULL, dat = NULL, Q = NULL, itemSelect = "GDI", MAXJ = 20, FIXED.LENGTH = TRUE, att.prior = NULL, post.initial = NULL, max.cut = 0.80,
                  NPS.args = list(gate = NULL, pseudo.prob = T, w.type = 1, seed = NULL), n.cores = 2, i.print = 250){

  # `%dopar%` <- foreach::`%dopar%` # define a local dopar

  #-------------------------
  # Gather data and objects
  #-------------------------

  est <- fit
  if(!is.null(est)){
    if(!is.null(est$extra$call)){ # package GDINA
      model <- ifelse(length(unique(est$model)) == 1, unique(est$model), "Combination")
      if(is.null(dat)){X <- est$options$dat} else {X <- dat}
      if(is.null(Q)){Q <- est$options$Q}
      LC.prob <- est$LC.prob
    } else { # package CDM
      model <- ifelse(length(unique(est$rule)) == 1, unique(est$rule), "Combination")
      if(is.null(dat)){X <- est$data} else {X <- dat}
      if(is.null(Q)){Q <- est$q.matrix}
      LC.prob <- matrix(NA, nrow = nrow(Q), ncol = 2^ncol(Q), dimnames = list(paste("Item", 1:nrow(Q)),apply(GDINA::attributepattern(ncol(Q)), 1, paste, collapse = "")))
      for(j in 1:nrow(Q)){
        kj <- which(Q[j,] == 1)
        prob.j.l <- est$probitem[est$probitem$itemno == j,]$prob
        tmp <- factor(apply(GDINA::attributepattern(ncol(Q))[,kj, drop = F], 1, paste, collapse = ""))
        tmp <- factor(tmp, levels = apply(unique(GDINA::attributepattern(length(kj))), 1, paste, collapse = ""))
        for(l in 1:length(prob.j.l)){LC.prob[j, which(as.numeric(tmp) == l)] <- prob.j.l[l]}
      }
    }
    LC.prob[LC.prob > 1] <- 1
    LC.prob[LC.prob < 0] <- 0
  } else { # nonparametric
    model <- paste0("NP_", NPS.args$gate)
    X <- dat
  }
  N <- nrow(X)
  J <- nrow(Q)
  K <- ncol(Q)
  L <- 2^K
  pattern <- GDINA::attributepattern(K)
  Lclass <- apply(GDINA::attributepattern(K), MARGIN = 1,FUN = function(x){paste(x,collapse = "")})
  if(is.null(att.prior)){att.prior <- rep(1/L, L)}
  if(is.null(post.initial)){post.initial <- matrix(data = rep(rep(1/L, L), N), nrow = N)}
  cl <- parallel::makeCluster(n.cores, type = "SOCK")
  doParallel::registerDoParallel(cl, cores = n.cores)

  #----------------------------
  # Warning end error messages
  #----------------------------

  if(!(itemSelect %in% c("GDI", "JSD", "MPWKL", "PWKL", "NPS", "random"))){stop("itemSelect must be one of the following: 'GDI', 'JSD', 'MPWKL', 'PWKL', 'NPS', 'random'")}
  if(itemSelect == "NPS"){
    if(is.null(dat)){stop("dat required when itemSelect == 'NPS'")}
    if(is.null(Q)){stop("Q required when itemSelect == 'NPS'")}
    if(is.null(NPS.args$gate)){
      stop("gate == 'AND' or gate == 'OR' required in NPS.args")
    } else {
      if(!(NPS.args$gate %in% c("AND", "OR"))){stop("gate == 'AND' or gate == 'OR' required in NPS.args")}
    }
    if(!(NPS.args$pseudo.prob %in% c(TRUE, FALSE))){stop("pseudo.prob == TRUE or pseudo.prob == FALSE required in NPS.args")}
    if(!FIXED.LENGTH & !NPS.args$pseudo.prob){warning("FIXED.LENGHT = FALSE is not available with itemSelect == 'NPS' if NPS.args$pseudo.prob == FALSE. NPS.args$pseudo.prob = TRUE applied instead")}
  } else {
    if(is.null(fit)){stop("fit required when itemSelect != 'NPS'")}
  }
  if(n.cores > parallel::detectCores()){stop("n.cores cannot be higher than the number of cores in the CPU")}

  #-----------------------
  # CD-CAT implementation
  #-----------------------

  out <- list()

  if(itemSelect != "NPS"){

    #-------------------
    # Parametric CD-CAT
    #-------------------

    out <- foreach(i = 1:N, .export = c("GDI.M", "H", "JSD.DICO.M", "PWKL.M", "MPWKL.M"), .inorder = T) %dopar% {
      try({
        if(i %% i.print == 0){cat(i)}

        # Start examinee i
        mlogPost_GDI <- post.initial[i,]
        est.GDI <- matrix(NA, nrow = 1, ncol = 7 + K, dimnames = list(1, c("ML", "nmodesML", "Lik", "MAP", "nmodesMAP", "Post", "EAP", paste("K", 1:K, sep = ""))))[-1,]
        item.log <- NULL

        # Item selection rule
        GDI <- switch(itemSelect,
                      "GDI" = GDI.M(LC.prob = LC.prob, mlogPost_GDI),
                      "JSD" = JSD.DICO.M(LC.prob = LC.prob, mlogPost_GDI),
                      "MPWKL" = MPWKL.M(LC.prob = LC.prob, mlogPost_GDI),
                      "PWKL" = PWKL.M(LC.prob = LC.prob, mlogPost_GDI, point.est = sample(1:L, size = 1)),
                      "random" = runif(J, 0, 1))
        names(GDI) <- 1:length(GDI)

        # Apply CD-CAT to examinee i
        jjselect <- 1
        jjcatGDI <- IP <- NULL
        for(jjselect in 1:MAXJ){
          if(is.null(jjcatGDI)){
            jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI)))
          } else {
            jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI[-c(jjcatGDI)])))
          }

          x.jj <- X[i, jjcatGDI[jjselect]]
          IP <- cbind(IP, x.jj)

          p.xi <- matrix(LC.prob[jjcatGDI,], ncol = L)
          p.xi.aj1_GDI <- matrix(0, nrow = L, ncol = length(IP))
          for(l in 1L:L){
            for(jj in 1L:length(IP)){
              res <- 1
              p.xi.aj1_GDI[l,jj] <- res * (p.xi[jj, l]^IP[jj]) * ((1 - p.xi[jj, l])^(1 - IP[jj]))
              res <- p.xi.aj1_GDI[l, jj]
            }
          }

          p.xi.aj_GDI <- (t(as.matrix(apply(p.xi.aj1_GDI, 1, prod))))
          mPost_GDI <- p.xi.aj_GDI * att.prior
          msumPost_GDI <- sum(mPost_GDI)
          mlogPost_GDI <- mPost_GDI / msumPost_GDI
          p.aj_GDI <- mlogPost_GDI %*% pattern

          # Item selection rule
          GDI <- switch(itemSelect,
                        "GDI" = GDI.M(LC.prob = LC.prob, mlogPost_GDI),
                        "JSD" = JSD.DICO.M(LC.prob = LC.prob, mlogPost_GDI),
                        "MPWKL" = MPWKL.M(LC.prob = LC.prob, mlogPost_GDI),
                        "PWKL" = PWKL.M(LC.prob = LC.prob, mlogPost_GDI, point.est = sample(1:L, size = 1)),
                        "NP" = stop("itemSelect == 'NP' is currently not available"),
                        "random" = runif(J, 0, 1))
          names(GDI) <- 1:length(GDI)

          est.GDI <- rbind(est.GDI,
                           c(Lclass[which.max(p.xi.aj_GDI)], length(which(p.xi.aj_GDI == max(p.xi.aj_GDI))), max(p.xi.aj_GDI),
                             Lclass[which.max(mlogPost_GDI)], length(which(mlogPost_GDI == max(mlogPost_GDI))), max(mlogPost_GDI),
                             paste(as.numeric(p.aj_GDI > .50), collapse = ""), p.aj_GDI))

          jjselect <- jjselect + 1
          if((jjselect) > MAXJ){break}
          if(!FIXED.LENGTH){if(max(mlogPost_GDI) >= max.cut){break}} # if FIXED.LENGTH == FALSE
        } # end CD-CAT for examinee i

        # Gather results for examinee i
        item.log <- c(item.log, jjcatGDI)
        out[[i]] <- list(est.cat.GDI = est.GDI, item.usage = item.log)
      })
    } # end parametric CD-CAT for all examinees

  } else {

    #-----------------------
    # Nonparametric CD-CAT
    #-----------------------

    att.pat <- matrix(NA, nrow = N, ncol = K) # Estimated attribute patterns
    exp.items <- list() # Items exponsed
    resp.pat <- list() # Examinees' response pattern
    if(NPS.args$gate == "AND"){gate <- rep("AND", J)}
    if(NPS.args$gate == "OR"){gate <- rep("OR", J)}
    if(!FIXED.LENGTH & !NPS.args$pseudo.prob){NPS.args$pseudo.prob <- TRUE}

    out <- foreach(i = 1:N, .packages = "NPCD", .export = c("pseudoP", "NPC.eta"), .inorder = T) %dopar% {
      try({
        if(i %% i.print == 0){cat(i)}
        if(!is.null(NPS.args$seed)){set.seed(NPS.args$seed + i)}

        if(NPS.args$pseudo.prob){
          out.i <- matrix(NA, ncol = 8 + K)[-1,]
          colnames(out.i) <- c("j", "qj", "yj", "alpha", "HD.alpha", "alpha2", "HD.alpha2", "HD.diff", paste0("pP.K", 1:K))
        } else {
          out.i <- matrix(NA, ncol = 8)[-1,]
          colnames(out.i) <- c("j", "qj", "yj", "alpha", "HD.alpha", "alpha2", "HD.alpha2", "HD.diff")
        }

        x <- X[i,] # fixed response vector for examinee i
        y <- c() # adaptive response vector for examinee i
        Qy <- c() # items presented to examinee i
        gatey <- c()

        # Step 0: Initialize item pool
        R0 <- as.data.frame(Q)
        gate0 <- gate

        # Step 1: Select first item with qj = e1, where e1 = 1 and all other entries = 0
        e1 <- c(1, rep(0, K - 1))
        hit <- which(sapply(apply(R0, 1, function(x) which(all(as.double(x) == e1))), sum) == 1)
        if(length(hit) > 1){hit <- sample(hit, 1)}
        R1 <- R0[-hit,]
        gate1 <- gate0[-hit]

        # Step 2: Administer item j (hit) to the examinee
        y <- c(y, x[hit])
        x <- x[-hit]
        Qy <- rbind(Qy, R0[hit,])
        gatey <- c(gatey, gate0[hit])

        if(NPS.args$pseudo.prob){
          out.i <- rbind(out.i, c(as.double(rownames(Qy)), paste(Qy[1,], collapse = ""), y, rep(NA, 5 + K)))
        } else {
          out.i <- rbind(out.i, c(as.double(rownames(Qy)), paste(Qy[1,], collapse = ""), y, rep(NA, 5)))
        }

        # Steps 3 to 5: Administer up to K items, based on the Q-optimal criterion for DINA gate (Xu, Wang, & Shang, 2016)
        for(k in 2:K){
          hit <- NULL

          # Step 3: Select new item and update R
          assign(paste0("e", k), c(rep(0, k - 1), 1, rep(0, K - k)))
          while(length(hit) == 0){
            qj.prev <- matrix(get(paste0("e", k)), ncol = K, byrow = T)
            for(w in 1:(k - 1)){qj.prev <- rbind(qj.prev, sample(c(0, 1), 1) * get(paste0("e", w)) * y[w])}
            qj <- colSums(qj.prev)
            qj[qj > 0] <- 1
            hit <- as.vector(which(sapply(apply(get(paste0("R", k - 1)), 1, function(x) which(all(as.double(x) == qj))), sum) == 1))
            if(length(hit) > 1){hit <- sample(hit, 1)} # hit is not equal to item number: it is related to R, not to Q
          }
          assign(paste0("R", k), get(paste0("R", k - 1))[-hit,])
          assign(paste0("gate", k), get(paste0("gate", k - 1))[-hit])

          # Step 4: Administer item j (hit) to the examinee
          y <- c(y, x[hit])
          x <- x[-hit]
          Qy <- rbind(Qy, get(paste0("R", k - 1))[hit,])
          gatey <- c(gatey, get(paste0("gate", k - 1))[hit])

          if(NPS.args$pseudo.prob){
            out.i <- rbind(out.i, c(as.double(rownames(Qy))[k], paste(Qy[k,], collapse = ""), y[k], rep(NA, 5 + K)))
          } else {
            out.i <- rbind(out.i, c(as.double(rownames(Qy))[k], paste(Qy[k,], collapse = ""), y[k], rep(NA, 5)))
          }

          # Step 5: k = k + 1 until k = K
        }

        # Step 6: Estimate the examinee's attribute profile (NPC method) and the second most plausible attribute profile
        # assign(paste0("a.", k), as.double(NPC(y, Qy, gatey)$HD))
        alphanp <- NPCD::AlphaNP(t(as.matrix(y)), Qy, gatey[1])
        # assign(paste0("a.", k), alphanp$alpha.est)
        # a2 <- names(sort(NPC(y, Qy, gatey)$HD.l)[m])
        pat.dist <- cbind(alphanp$pattern, alphanp$loss.matrix)
        if(K > 1){random.order <- sample(1:K, K, F)}
        for(k in random.order){pat.dist <- pat.dist[order(pat.dist[,k]),]}
        pat.dist <- pat.dist[order(pat.dist[, K + 1]),]
        assign(paste0("a.", K), pat.dist[order(pat.dist[,K + 1]),][1, 1:K])
        m <- 2
        assign(paste0("a2.", K), pat.dist[order(pat.dist[,K + 1]),][m, 1:K])

        out.i[K, 4:8] <- c(paste(get(paste0("a.", K)), collapse = ""), pat.dist[order(pat.dist[,K + 1]),][1, K + 1],
                           paste(get(paste0("a2.", K)), collapse = ""), pat.dist[order(pat.dist[,K + 1]),][2, K + 1],
                           pat.dist[order(pat.dist[,K + 1]),][2, K + 1] - pat.dist[order(pat.dist[,K + 1]),][1, K + 1])

        if(NPS.args$pseudo.prob){
          pP.k <- as.vector(pseudoP(alphanp, inv.type = "min", w.type = NPS.args$w.type))
          out.i[K, 9:(8 + K)] <- round(pP.k, 5)
        }

        # Step 7: Select item to administer that can discriminate a.k and a2.k
        t <- K + 1
        while(t <= MAXJ){
          hit <- NULL
          while(is.null(hit)){
            items <- sample(1:nrow(get(paste0("R", t - 1))))
            for(j in items){
              qj <- as.double(get(paste0("R", t - 1))[j,])
              eta.a <- NPC.eta(qj, get(paste0("a.", t - 1)), get(paste0("gate", t - 1))[j])
              eta.a2 <- NPC.eta(qj, get(paste0("a2.", t - 1)), get(paste0("gate", t - 1))[j])
              if(eta.a != eta.a2){

                # Step 10: Administer item and update y and R
                hit <- as.vector(which(sapply(apply(get(paste0("R", t - 1)), 1, function(x) which(all(as.double(x) == qj))), sum) == 1))
                if(length(hit) > 1){hit <- sample(hit, 1)} # hit is not equal to item number: it is related to R, not to Q
                assign(paste0("R", t), get(paste0("R", t - 1))[-hit,])
                assign(paste0("gate", t), get(paste0("gate", t - 1))[-hit])
                y <- c(y, x[hit])
                x <- x[-hit]
                Qy <- rbind(Qy, get(paste0("R", t - 1))[hit,])
                gatey <- c(gatey, get(paste0("gate", t - 1))[hit])
                m <- 2
                break
              } else {
                next
              }
            }

            # Step 8 and 9: Find an item that can discriminate a.k and a2.k
            if(is.null(hit)){
              m <- m + 1
              # a2 <- names(sort(NPC(y, Qy, gatey)$HD.l)[m])
              alphanp <- NPCD::AlphaNP(t(as.matrix(y)), Qy, gatey[1])
              pat.dist <- cbind(alphanp$pattern, alphanp$loss.matrix)
              if(K > 1){random.order <- sample(1:K, K, F)}
              for(k in random.order){pat.dist <- pat.dist[order(pat.dist[,k]),]}
              pat.dist <- pat.dist[order(pat.dist[, K + 1]),]
              assign(paste0("a2.", t - 1), pat.dist[order(pat.dist[,K + 1]),][m, 1:K])
            }
          }

          # Step 11: Estimate the examinee's attribute profile (NPC method) and the second most plausible attribute profile
          # assign(paste0("a.", k), as.double(NPC(y, Qy, gatey)$HD))
          alphanp <- NPCD::AlphaNP(t(as.matrix(y)), Qy, gatey[1])
          # assign(paste0("a.", t), alphanp$alpha.est)
          # a2 <- names(sort(NPC(y, Qy, gatey)$HD.l)[m])
          pat.dist <- cbind(alphanp$pattern, alphanp$loss.matrix)
          if(K > 1){random.order <- sample(1:K, K, F)}
          for(k in random.order){pat.dist <- pat.dist[order(pat.dist[,k]),]}
          pat.dist <- pat.dist[order(pat.dist[, K + 1]),]
          assign(paste0("a.", t), pat.dist[order(pat.dist[,K + 1]),][1, 1:K])
          m <- 2
          assign(paste0("a2.", t), pat.dist[order(pat.dist[,K + 1]),][m, 1:K])

          if(!NPS.args$pseudo.prob){
            out.i <- rbind(out.i, c(as.double(rownames(Qy))[t], paste(Qy[t,], collapse = ""), y[t],
                                    paste(get(paste0("a.", t)), collapse = ""), pat.dist[order(pat.dist[,K + 1]),][1, K + 1],
                                    paste(get(paste0("a2.", t)), collapse = ""), pat.dist[order(pat.dist[,K + 1]),][2, K + 1],
                                    pat.dist[order(pat.dist[,K + 1]),][2, K + 1] - pat.dist[order(pat.dist[,K + 1]),][1, K + 1]))
          } else {
            if(NPS.args$pseudo.prob){
              pP.k <- as.vector(pseudoP(alphanp, inv.type = "min", w.type = NPS.args$w.type))
              out.i <- rbind(out.i, c(as.double(rownames(Qy))[t], paste(Qy[t,], collapse = ""), y[t],
                                      paste(get(paste0("a.", t)), collapse = ""), pat.dist[order(pat.dist[,K + 1]),][1, K + 1],
                                      paste(get(paste0("a2.", t)), collapse = ""), pat.dist[order(pat.dist[,K + 1]),][2, K + 1],
                                      pat.dist[order(pat.dist[,K + 1]),][2, K + 1] - pat.dist[order(pat.dist[,K + 1]),][1, K + 1],
                                      round(pP.k, 5)))
            }
          }

          if(NPS.args$pseudo.prob){ # if FIXED.LENGTH == FALSE
            if(!FIXED.LENGTH){
              cut.pP.k <- pP.k
              cut.pP.k[cut.pP.k < 0.5] <- 1 - cut.pP.k[cut.pP.k < 0.5]
              if(all(cut.pP.k > max.cut)){break}
            }
          }

          t <- t + 1
        }
        att.pat[i,] <- get(paste0("a.", t - 1))
        exp.items[[i]] <- as.double(rownames(Qy))
        out.i <- as.data.frame(out.i)
        out.i$j <- as.double(as.character(out.i$j))
        out.i$qj <- as.character(out.i$qj)
        out.i$yj <- as.double(as.character(out.i$yj))
        out.i$alpha <- as.character(out.i$alpha)
        out.i$HD.alpha <- as.double(as.character(out.i$HD.alpha))
        out.i$alpha2 <- as.character(out.i$alpha2)
        out.i$HD.alpha2 <- as.double(as.character(out.i$HD.alpha2))
        out.i$HD.diff <- as.double(as.character(out.i$HD.diff))
        if(NPS.args$pseudo.prob){for(k in 1:K){out.i[,8 + k] <- as.double(as.character(out.i[,8 + k]))}}
        out[[i]] <- list(est.cat.NPS = out.i, item.usage = exp.items[[i]])
      })
    } # end nonparametric CD-CAT for all examinees
  } # end CD-CAT

  parallel::stopCluster(cl)
  res <- list()
  res$est <- out
  res$specifications <- list("fit" = fit, "dat" = dat, "Q" = Q, "model" = model,
                             "itemSelect" = itemSelect, "MAXJ" = MAXJ, "FIXED.LENGTH" = FIXED.LENGTH,
                             "att.prior" = att.prior, "post.initial" = post.initial, "max.cut" = max.cut,
                             "NPS.args" = NPS.args, "n.cores" = n.cores, "i.print" = i.print)
  class(res) <- "cdcat"
  return(res)
}
