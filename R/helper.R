proggresive.f <- function(jjselect, MAXJ, b, GDI){
  # used in cdcat function
  # performs progressive method on itemSelect (Revuelta & Ponsoda, 1998; Barrada, Olea, Ponsoda, & Abad, 2008)
  if(jjselect == 1){
    W_jjGDI <- 0 
  }else{
    W_jjGDI <- sum((c(2:jjselect)-1)^b)/ sum((c(2:MAXJ)-1)^b)
  }
  R_jjGDI <- runif(length(GDI),0,max(GDI)) # length(GDI) por j 22/06
  S.GDI.jj <- ((1-rep(W_jjGDI,length(R_jjGDI)))*R_jjGDI + rep(W_jjGDI, length(R_jjGDI))*GDI)
  return(S.GDI.jj)
}

overlap.rate <- function(CAT.length, exposure.rates){
  # used in cdcat.getdata function
  Q <- CAT.length
  exp <- exposure.rates
  n <- length(exp)
  # Q: test length, n: Item bank size, exp: Item exposure rates
  res <- (n/Q)*var(exp) + Q/n # Barrada et al. (2009) Psicothema
  return(res)
}

cdcat.getdata <- function(cdcat.obj, alpha = NULL, plots = TRUE){
  # used in cdcat.summary function
  # Gather information
  if(!inherits(cdcat.obj, "cdcat")){stop("cdcat class required.")}
  dat <- cdcat.obj$specifications$dat
  if(is.null(dat)){dat <- cdcat.obj$specifications$fit$options$dat}
  Q <- cdcat.obj$specifications$Q
  if(is.null(Q)){Q <- cdcat.obj$specifications$fit$options$Q}
  N <- nrow(dat)
  J <- nrow(Q)
  K <- ncol(Q)
  MAXJ <- cdcat.obj$specifications$MAXJ
  FIXED.LENGTH <- cdcat.obj$specifications$FIXED.LENGTH
  itemSelect <- cdcat.obj$specifications$itemSelect
  
  # Gather attribute profiles and posterior probabilities
  alpha.est <- alpha.pp <- matrix(NA, nrow = N, ncol = K, dimnames = list(NULL, colnames(Q)))
  item.exp <- list()
  for(i in 1:N){
    tmp <- cdcat.obj$est[[i]]$est.cat
    if(itemSelect %in% c("NPS", "GNPS")){
      alpha.est[i,] <- as.numeric(unlist(strsplit(tmp$alpha[nrow(tmp)], "", "")))
      if(cdcat.obj$specifications$NP.args$PPP){
        alpha.pp[i,] <- as.numeric(tmp[nrow(tmp), (ncol(tmp) - K + 1):ncol(tmp)])
      }
    } else {
      alpha.est[i,] <- as.numeric(unlist(strsplit(tmp$MAP[nrow(tmp)], "", "")))
      alpha.pp[i,] <- as.numeric(tmp[nrow(tmp), (ncol(tmp) - K + 1):ncol(tmp)])
    }
    item.exp[[i]] <- cdcat.obj$est[[i]]$item.usage
  }
  
  # Attribute mastery
  att.prevalence <- data.frame(Attribute = names(colMeans(alpha.est)),
                               Prevalence = as.numeric(colMeans(alpha.est)))
  att.prevalence$Attribute <- factor(att.prevalence$Attribute, levels = att.prevalence$Attribute)
  if(plots){
    att.prevalence.plot <- ggplot2::ggplot(att.prevalence, ggplot2::aes(x = Attribute, y = Prevalence)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      ggplot2::theme_bw(base_size = 15) +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
  att.profiles <- data.frame(Profile = names(table(apply(alpha.est, 1, paste, collapse = ""))),
                             Frequency = as.numeric(table(apply(alpha.est, 1, paste, collapse = ""))))
  att.profiles$Profile <- factor(att.profiles$Profile, levels =   att.profiles[order(att.profiles$Frequency, decreasing = TRUE),]$Profile)
  if(plots){
    att.profiles.plot <- ggplot2::ggplot(att.profiles, ggplot2::aes(x = Profile, y = Frequency)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_y_continuous() +
      ggplot2::theme_bw(base_size = 15) +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))
  }
  
  # Attribute probability
  if(!(itemSelect %in% c("NPS", "GNPS") & !cdcat.obj$specifications$NP.args$PPP)){
    att.prob <- data.frame(Attribute = rep(colnames(alpha.pp), each = N),
                           Probability = as.numeric(alpha.pp))
    att.prob$Attribute <- factor(att.prob$Attribute, levels = unique(att.prob$Attribute))
    if(plots){
      att.prob.plot <- ggplot2::ggplot(att.prob, ggplot2::aes(x = Attribute, y = Probability, fill = Probability)) +
        ggplot2::geom_violin(col = "transparent", fill = "black", alpha = 0.05) +
        ggplot2::geom_point(size = 2.2, pch = 21, col = "black", alpha = 0.3) +
        ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
        ggplot2::scale_fill_gradient(low = "firebrick3", high = "seagreen3") +
        ggplot2::theme_bw(base_size = 15) +
        ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::guides(fill = "none")
    }
  }
  
  # Item exposition
  test.length <- data.frame(Length = as.numeric(names(table(sapply(item.exp, length)))),
                            Frequency = as.numeric(table(sapply(item.exp, length))))
  if(plots){
    test.length.plot <- ggplot2::ggplot(test.length, ggplot2::aes(x = Length, y = Frequency)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_x_continuous("Test Length", limits = c(0, MAXJ + 1), breaks = seq(0, MAXJ, 1)) +
      ggplot2::theme_bw(base_size = 15) +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank())
  }
  
  tmp <- table(unlist(lapply(cdcat.obj$est, function(x) x[[2]])))
  df.item.exp <- data.frame(Item = 1:J, Exp = 0)
  df.item.exp[as.numeric(names(tmp)), 2] <- tmp/N
  if(plots){
    df.item.exp.plot <- ggplot2::ggplot(df.item.exp, ggplot2::aes(x = Item, y = Exp)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_x_continuous("Item", limits = c(0, max(df.item.exp$Item) + 1), breaks = seq(0, max(df.item.exp$Item), floor(max(df.item.exp$Item) / 10))) +
      ggplot2::scale_y_continuous("Exposure Rate", limits = c(0, 1)) +
      ggplot2::theme_bw(base_size = 15) +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank())
  }
  item.exp.stats <- c(summary(df.item.exp$Exp), `overlap rate` = overlap.rate(weighted.mean(test.length$Length, test.length$Frequency), df.item.exp$Exp))
  
  # Classification accuracy
  if(!is.null(alpha)){
    est <- lapply(cdcat.obj$est, "[[", 1)
    df.recovery <- data.frame(ItemPos = 1:MAXJ, PCV = NA, PCA = NA)
    for(l in 1:MAXJ){
      if(itemSelect %in% c("NPS", "GNPS")){
        pre.tmp <- lapply(lapply(lapply(est, function(x) x[l, "alpha"]), as.character), FUN = strsplit, split = "")
      } else {
        pre.tmp <- lapply(lapply(lapply(est, function(x) x[l, "MAP"]), as.character), FUN = strsplit, split = "")
      }
      pre.tmp[which(sapply(pre.tmp, function(x) any(is.na(x))))] <- as.list(as.data.frame(matrix(NA, nrow = K, ncol = length(which(sapply(pre.tmp, function(x) any(is.na(x))))))))
      # for(nl in which(sapply(pre.tmp, function(x) any(is.na(x))))){
      #   pre.tmp[[nl]] <- rep(NA, K)
      # }
      tmp <- matrix(data = as.numeric(unlist(pre.tmp)), ncol = K, nrow = N, byrow = TRUE)
      tmp[-which(complete.cases(tmp)),] <- alpha.est[-which(complete.cases(tmp)),]
      df.recovery$PCV[l] <- mean(rowSums(tmp == alpha) == K)
      df.recovery$PCA[l] <- mean(tmp == alpha)
    }
    PCA <- c(sapply(1:K, function(k) mean(alpha.est[,k] == alpha[,k])), GDINA::ClassRate(alpha.est, alpha)$PCA)
    names(PCA) <- c(paste0("K", 1:K), "mean.PCA")
    PCV <- GDINA::ClassRate(alpha.est, alpha)$PCV
    names(PCV) <- paste0(1:K, "/", K)
    if(plots){
      plotPCAm <- ggplot2::ggplot(df.recovery, ggplot2::aes(x = ItemPos, y = PCA)) +
        ggplot2::geom_line() +
        ggplot2::geom_point(size = 2) +
        ggplot2::scale_y_continuous("Attribute Recovery", limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
        ggplot2::theme_bw(base_size = 15) +
        ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(), panel.grid.minor.y = ggplot2::element_blank())
      plotPCV <- ggplot2::ggplot(df.recovery, ggplot2::aes(x = ItemPos, y = PCV)) +
        ggplot2::geom_line() +
        ggplot2::geom_point(size = 2) +
        ggplot2::scale_y_continuous("Pattern Recovery", limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
        ggplot2::theme_bw(base_size = 15) +
        ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(), panel.grid.minor.y = ggplot2::element_blank())
      if(FIXED.LENGTH){
        plotPCAm <- plotPCAm + ggplot2::scale_x_continuous("Until Item Position", breaks = seq(0, MAXJ, 1))
        plotPCV <- plotPCV + ggplot2::scale_x_continuous("Until Item Position", breaks = seq(0, MAXJ, 1))
      } else {
        plotPCAm <- plotPCAm + ggplot2::scale_x_continuous("Until (Accumulated) Item Position", breaks = seq(0, MAXJ, 1))
        plotPCV <- plotPCV + ggplot2::scale_x_continuous("Until (Accumulated) Item Position", breaks = seq(0, MAXJ, 1))
      }
    }
  }
  
  if(itemSelect == "GNPS" | !cdcat.obj$specifications$NP.args$PPP){
    alpha.pp <- att.prob <- att.prob.plot <- NULL
  }
  if(is.null(alpha)){
    df.recovery <- PCA <- plotPCAm <- PCV <- plotPCV <- NULL
  }
  if(plots){
    alpha.estimates <- list(est = alpha.est, pp = alpha.pp, 
                            prevalence = att.prevalence, prevalence.plot = att.prevalence.plot,
                            post = att.prob, post.plot = att.prob.plot,
                            profiles = att.profiles, profiles.plot = att.profiles.plot)
    item.exposure <- list(stats = item.exp.stats,
                          exp = df.item.exp, exp.plot = df.item.exp.plot,
                          length = test.length, length.plot = test.length.plot)
    alpha.recovery <- list(recovery = df.recovery,
                           PCA = PCA, PCA.plot = plotPCAm,
                           PCV = PCV, PCV.plot = plotPCV)
  } else {
    alpha.estimates <- list(est = alpha.est, pp = alpha.pp, 
                            prevalence = att.prevalence, post = att.prob, profiles = att.profiles)
    item.exposure <- list(stats = item.exp.stats,
                          exp = df.item.exp, length = test.length)
    alpha.recovery <- list(recovery = df.recovery,
                           PCA = PCA, PCV = PCV)
  }
  if(is.null(alpha)){
    res <- list(alpha.estimates = alpha.estimates, item.exposure = item.exposure, specifications = cdcat.obj$specifications)
  } else {
    res <- list(alpha.estimates = alpha.estimates, item.exposure = item.exposure, alpha.recovery = alpha.recovery, specifications = cdcat.obj$specifications)
  }
  class(res) <- "cdcat.summary"
  return(res)
}
print.cdcat.summary <- function(x){
  if(length(x$specifications) == 1){
    models <- rownames(x$recovery$PCV)
    mss <- "========================================================="
    if(is.null(x$recovery$PCVcomp)){
      mss <- paste0(mss, "\n", "Variable length CD-CAT comparison of models: ", paste(models, collapse = ", "))
    } else {
      mss <- paste0(mss, "\n", "Fixed length CD-CAT comparison of models: ", paste(models, collapse = ", "))
    }
    mss <- paste0(mss, "\n", "\n",
                  "\t", paste(models, collapse = "\t"), "\n",
                  "PCA", "\t", paste(x$recovery$PCA[,ncol(x$recovery$PCA)], collapse = "\t"), "\n",
                  "PCV", "\t", paste(x$recovery$PCV[,ncol(x$recovery$PCV)], collapse = "\t"), "\n",
                  "Overlap", "\t", paste(round(x$item.exposure$stats[,ncol(x$item.exposure$stats)], 3), collapse = "\t"), "\n")
    mss <- paste0(mss, "=========================================================")
  } else {
    itemSelect <- x$specifications$itemSelect
    FIXED.LENGTH <- x$specifications$FIXED.LENGTH
    MAXJ <- x$specifications$MAXJ
    att.prevalence <- round(colMeans(x$alpha.estimates$est), 3)
    MEANJ <- weighted.mean(x$item.exposure$length[,1], x$item.exposure$length[,2])
    overlap.rate <- x$item.exposure$stats["overlap rate"]
    mss <- "========================================================="
    if(FIXED.LENGTH){
      mss <- paste0(mss, "\n", "Fixed length CD-CAT")
    } else {
      mss <- paste0(mss, "\n", "Variable length CD-CAT")
    }
    mss <- paste0(mss, " using the ", itemSelect, " item selection rule", "\n", "\n")
    mss <- paste0(mss, 
                  "Maximum number of items = ", MAXJ, "\n", 
                  "Average number of items = ", MEANJ, "\n",
                  "Overlap rate = ", round(overlap.rate, 3), "\n","\n",
                  "Attribute prevalence:", "\n")
    for(k in 1:length(att.prevalence)){mss <- paste0(mss, names(att.prevalence)[k], " = ", round(as.numeric(att.prevalence[k]), 3), "\n")}
    if(!is.null(x$alpha.recovery)){
      mss <- paste0(mss, "\n",
                    "Attribute class. accuracy (PCA) = ", round(as.numeric(x$alpha.recovery$PCA["mean.PCA"]), 3), "\n",
                    "Pattern class. accuracy (PCV) = ", round(as.numeric(x$alpha.recovery$PCV[length(x$alpha.recovery$PCV)]), 3), "\n")
    }
    mss <- paste0(mss, "=========================================================")
  }
  cat(mss)
}

genQ <- function(J, K, n.id = 2, qkProp, PropWithId = T, min.jk = 1, max.kcor = 0.3, seed = NULL){
  # used in gen.itembank function
  # J: number of items
  # K: number of attributes
  # qkProp: proportion of non-identity matrix items with k attributes (e.g., qkProb = c(0.1, 0.45, 0.45), 0.1 of items will have 1 attribute, 0.45 will have 2 attributes, etc.)
  # n.id: number of identity matrices inside Q-matrix
  # min.jk: minimum number of non-identity matrix items measuring each attribute
  # max.kcor: maximum permitted attribute correlation
  # seed: seed

  if(!is.null(seed)){set.seed(seed)}
  if(length(qkProp) > K){stop("length(qkProp) must be <= K")}
  if (!isTRUE(all.equal(1, sum(qkProp), tolerance = .Machine$double.eps^0.25))) {stop("sum(qkProp) must be equal to 1")}
  pat <- GDINA::attributepattern(K)
  for(k in 1:length(qkProp)){assign(paste0("pat", k), rbind(pat[which(rowSums(pat) == k),]))}
  if((K * n.id) > J){stop("Identity matrix/matrices has/have more rows than the whole Q-matrix ((K * nI) > J)")}
  I <- matrix(rep(diag(1, K), n.id), ncol = K, byrow = T)
  Q <- c()

  if(!PropWithId){
    Jq <- J - K * n.id
    Jk <- c()
    for(k in 1:length(qkProp)){Jk <- c(Jk, floor(Jq * qkProp[k]))}
    while(sum(Jk) < Jq){
      tmp <- sample(1:length(Jk), 1)
      Jk[tmp] <- Jk[tmp] + 1
    }
    for(k in 1:length(qkProp)){assign(paste0("replace", k), ifelse(Jk[k] > nrow(get(paste0("pat", k))), T, F))}
    Q <- rbind(I)
    for(k in 1:length(qkProp)){Q <- rbind(Q, rbind(get(paste0("pat", k))[sample(nrow(get(paste0("pat", k))), size = Jk[k], replace = get(paste0("replace", k))),]))}
    i <- 0
    while(any(colSums(Q) < min.jk) | any(cor(Q)[lower.tri(cor(Q))] > max.kcor)){
      Q <- rbind(I)
      for(k in 1:length(qkProp)){Q <- rbind(Q, rbind(get(paste0("pat", k))[sample(nrow(get(paste0("pat", k))), size = Jk[k], replace = get(paste0("replace", k))),]))}
      i <- i + 1
      if(i == 500){stop("In gen.itembank, no convergence in Q-matrix generation. Try a lower value for argument gen.Q$minJ.K or a higher value for argument gen.Q$max.Kcor")}
    }
  } else {
    Jk <- c()
    for(k in 1:length(qkProp)){Jk <- c(Jk, floor(J * qkProp[k]))}
    while(sum(Jk) < J){
      tmp <- sample(which(Jk != 0), 1)
      Jk[tmp] <- Jk[tmp] + 1
    }
    if(min.jk > (sum(Jk * 1:length(Jk)) / (J*K)) * J){cat("Warning in genQ: min.jk may be to high according to the other arguments. Consider using a lower value for min.jk")}
    for(k in 1:length(qkProp)){assign(paste0("replace", k), ifelse(Jk[k] > nrow(get(paste0("pat", k))), T, F))}
    Q <- rbind(I)
    Jk[1] <- Jk[1] - nrow(Q)
    for(k in 1:length(qkProp)){Q <- rbind(Q, rbind(get(paste0("pat", k))[sample(nrow(get(paste0("pat", k))), size = Jk[k], replace = get(paste0("replace", k))),]))}
    i <- 0
    while(any(colSums(Q) < min.jk) | any(cor(Q)[lower.tri(cor(Q))] > max.kcor)){
      Q <- rbind(I)
      for(k in 1:length(qkProp)){Q <- rbind(Q, rbind(get(paste0("pat", k))[sample(nrow(get(paste0("pat", k))), size = Jk[k], replace = get(paste0("replace", k))),]))}
      i <- i + 1
      if(i == 500){stop("In gen.itembank, no convergence in Q-matrix generation. Try a lower value for argument gen.Q$minJ.K or a higher value for argument gen.Q$max.Kcor")}
    }
  }

  Jk[1] <- Jk[1] + K * n.id
  Jk <- matrix(c(Jk, K*n.id), nrow = 1, dimnames = list("n.items", c(paste0("q", 1:length(qkProp)), "id")))
  att.cor <- cor(Q)
  options <- list(J = J, K = K, n.id = n.id,
                  qkProp = matrix(qkProp, nrow = 1, dimnames = list("qkProp", paste0("q", 1:length(qkProp)))),
                  PropWithId = PropWithId, min.jk = min.jk, max.kcor = max.kcor, seed = seed)

  return(list(Q = Q, Jk = Jk, att.cor = att.cor, options = options))
}

gen.ACDMparam <- function(P0, P1, k.j, min.delta = 0, seed = NULL){
  # used in gen.itembank function
  if((P1 - P0) / k.j < min.delta){
    warning("min.delta not achievable with P0 and P1. min.delta has been forced to be lower.")
    min.delta <- (P1 - P0) / k.j
  }
  if(!is.null(seed)){set.seed(seed)}
  if(k.j == 1){
    catprob.parm <- c("P(0)" = P0, "P(1)" = P1)
    delta.parm <- c(d0 = P0, d1 = P1 - P0)
    return(list(catprob.parm = catprob.parm, delta.parm = delta.parm))
  } else {
    delta.parm <- -1
    while(any(delta.parm < min.delta)){
      tmp <- runif(k.j, round(min.delta, 10), round(P1 - P0 - (min.delta * (k.j - 1)), 10))
      delta.parm <- tmp / sum(tmp) * (P1 - P0)
    }
    delta.parm <- c(P0, delta.parm)
    names(delta.parm) <- paste0("d", 0:k.j)
    lc <- GDINA::attributepattern(k.j)
    catprob.parm <- rep(P0, nrow(lc))
    names(catprob.parm) <- paste0("P(", apply(lc, 1, paste, collapse = ""), ")")
    for(k in 1:k.j){catprob.parm[lc[,k] == 1] <- catprob.parm[lc[,k] == 1] + delta.parm[k + 1]}
    return(list(catprob.parm = catprob.parm, delta.parm = delta.parm))
  }
}
gen.GDINAparam <- function(P0, P1, k.j, min.delta = 0, seed = NULL){
  # used in gen.itembank function
  if(!is.null(seed)){set.seed(seed)}
  if(k.j == 1){
    catprob.parm <- c("P(0)" = P0, "P(1)" = P1)
    delta.parm <- c(d0 = P0, d1 = P1 - P0)
    return(list(catprob.parm = catprob.parm, delta.parm = delta.parm))
  } else {
    lc <- GDINA::attributepattern(k.j)
    pos.k <- apply(lc, 2, function(x) which(x == 1))
    delta.k <- rep(-1, k.j)
    while(any(delta.k < min.delta)){
      catprob.parm <- delta.parm <- P0
      for(l in 2:(nrow(lc) - 1)){
        k <- which(pos.k == l, arr.ind = T)[,2]
        if(length(k) == 1){
          catprob.parm <- c(catprob.parm, runif(1, P0, P1))
          delta.parm <- c(delta.parm, catprob.parm[l] - P0)
        } else {
          pos <- pos.k[,k][pos.k[,k] < l]
          pos <- pos[rowSums(lc[pos,]) < length(k)]
          catprob.parm <- c(catprob.parm, runif(1, max(catprob.parm[pos]), P1))
          delta.parm <- c(delta.parm, catprob.parm[l] - P0 - sum(delta.parm[pos]))
        }
      }
      catprob.parm <- c(catprob.parm, P1)
      delta.parm <- c(delta.parm, P1 - sum(delta.parm))
      delta.k <- delta.parm[1 + 1:k.j]
    }
    return(list(catprob.parm = catprob.parm, delta.parm = delta.parm))
  }
}

# Item selection procedures for cdcat function
GDI.M <- function(LC.prob, prior){
  GDI.est <- c(0)
  for (jj in 1:nrow(LC.prob)) {
    MPj <- sum(LC.prob[jj,]*prior)
    GDI.est[jj] <- sum(prior*(LC.prob[jj,]-MPj)^2)    }
  return(GDI.est)
}
H <- function(v){
  v <- v[v > 0]
  return(sum(-v * log2(v)))
}
JSD.DICO.M <- function(LC.prob, prior){
  # prob.jlist: LC.prob for correct/incorrect and item
  # prior: attribute joint distribution that is used as weight
  # date: 06/29/2017
  JSD.est <- c(0)
  w <- prior
  for (jj in 1:nrow(LC.prob)) {
    p <-  rbind(LC.prob[jj, ], 1 - LC.prob[jj, ])
    w <- matrix(data = w, ncol = 1)
    JSD.est[jj] <- H(p %*% w) - apply(p, 2, H) %*% w
  }
  return(JSD.est)
}
PWKL.M <- function(LC.prob, prior, point.est){
  # LC.prob: P(alpha=1) for each item
  # prior:   attribute joint distribution
  # date: 29/06/2017
  PWKL.est <- c(0)
  for (jj in 1:nrow(LC.prob)) {
    res <- sum(log(LC.prob[jj, point.est]/LC.prob[jj, ])*LC.prob[jj, point.est]*prior)
    PWKL.est[jj] <- res + sum(log((1-LC.prob[jj, point.est])/(1-LC.prob[jj, ]))*(1-LC.prob[jj, point.est])*prior)
  }
  return(PWKL.est)
}
MPWKL.M <- function(LC.prob, prior){
  # MPWKL: Modified Posterior Weighted KL from Kaplan et al. (2015)
  # LC.prob: P(alpha=1) for each item
  # prior:   attribute joint distribution
  # date: 01/03/2017
  MPWKL.est0 <- c(0)
  MPWKL.est <- c(0)
  for (jj in 1:nrow(LC.prob)) {
    for (point.est in 1:ncol(LC.prob)) {
      res <- sum(log(LC.prob[jj, point.est]/LC.prob[jj, ])*LC.prob[jj, point.est]*prior)
      MPWKL.est0[point.est] <- res + sum(log((1-LC.prob[jj, point.est])/(1-LC.prob[jj, ]))*(1-LC.prob[jj, point.est])*prior)
    } # end point est.
    MPWKL.est[jj] <- sum(MPWKL.est0 * prior)
  } # end jj
  return(MPWKL.est)
}
# Nonparametric CD-CAT
NPC.eta <- function(q, l, gate){
  eta_k <- c()
  K <- length(q)
  if(gate == "AND"){for(k in 1:K){eta_k <- c(eta_k, l[k]^q[k]); eta <- prod(eta_k)}}
  if(gate == "OR"){for(k in 1:K){eta_k <- c(eta_k, (1 - l[k])^q[k]); eta <- 1 - prod(eta_k)}}
  return(eta)
}
PPP <- function(dist, pattern, w, pp = NULL, alpha.est = NULL){
  if(nrow(dist) != nrow(pattern)){stop("dist and pattern must have the same number of rows")}
  N <- ncol(dist)
  K <- ncol(pattern)
  if(is.null(N)){N <- 1}
  ppp <- list()
  RMSE <- CA <- rep(NA, length(w)); names(RMSE) <- names(CA) <- w
  for(x in 1:length(w)){
    W <- t(apply(dist, 2, function(i) (1/w[x])^(i - max(i))))
    ppp[[x]] <- sapply(1:K, function(k) sapply(1:N, function(i) weighted.mean(pattern[, k], W[i,])))
    if(!is.null(RMSE)){RMSE[x] <- sqrt(mean((pp - ppp[[x]])^2))}
    if(!is.null(alpha.est)){CA[x] <- ClassRate(round(ppp[[x]]), alpha.est)$PCA}
  }
  names(ppp) <- w
  if(!is.null(pp)){
    sel.ppp <- ppp[[which.min(RMSE)[1]]]
    res <- list(sel.PPP = sel.ppp, sel.w = as.numeric(names(which.min(RMSE))))
  } else if(!is.null(alpha.est)){
    sel.ppp <- ppp[[which.max(CA)[1]]]
    res <- list(sel.PPP = sel.ppp, sel.w = as.numeric(names(which.max(CA))))
  } else {
    res <- list()
  }
  if(length(w) == 1){
    res$PPP <- ppp[[1]]
  } else {
    res$PPP <- ppp
  }
  if(!is.null(pp)){
    res$RMSE.w <- RMSE
  } else if(!is.null(alpha.est)){
    res$CA.w <- CA
  }
  return(res)
}
