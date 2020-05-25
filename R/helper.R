overlap.rate <- function(CAT.length, exposure.rates){
  # used in cdcat.getdata function
  Q <- CAT.length
  exp <- exposure.rates
  n <- length(exp)
  # Q: test length, n: Item bank size, exp: Item exposure rates
  res <- (n/Q)*var(exp) + Q/n # Barrada et al. (2009) Psicothema
  return(res)
}

cdcat.getdata <- function(cdcat.obj, alpha){
  # used in cdcat.summary function
  MAXJ <- cdcat.obj$specifications$MAXJ
  if(cdcat.obj$specifications$itemSelect != "NPS"){
    if(cdcat.obj$specifications$FIXED.LENGTH == TRUE){
      est <- lapply(cdcat.obj$est, '[[', 1)
      N <- length(est)
      K <- ncol(cdcat.obj$specifications$Q)
      est.MAP.jj <- list()
      PCV.jj <- as.data.frame(matrix(0, nrow = MAXJ, ncol = 2))
      PCAm.jj <- as.data.frame(matrix(0, nrow = MAXJ, ncol = 2))
      colnames(PCV.jj) <- c("item.position", "pattern.recovery")
      colnames(PCAm.jj) <- c("item.position", "attribute.recovery")
      for (jj in 1:MAXJ) {
        est.MAP.jj[[jj]] <- matrix(data = as.numeric(unlist(lapply(lapply(lapply(est, function(x) x[jj, "MAP"]), as.character), FUN = strsplit, split = ""))),
                                   ncol = K, nrow = N, byrow = TRUE)
        PCV.jj[jj, ] <- c(jj, mean(rowSums(est.MAP.jj[[jj]] == alpha) == K))
        PCAm.jj[jj, ] <- c(jj, mean(est.MAP.jj[[jj]] == alpha))
      }
      recovery <- list()
      for (kk in 1:K) {
        recovery$PCV[kk] <- mean(rowSums(est.MAP.jj[[jj]] == alpha) >= kk)}
      recovery$PCA <- c(colMeans(est.MAP.jj[[jj]] == alpha), mean(est.MAP.jj[[jj]] == alpha))
      names(recovery$PCV) <- paste(1:K, sep = "/", K)
      names(recovery$PCA) <- c(paste("K", 1:K, sep = ""), "mean.PCA")
      recovery$dataplot <- cbind(PCV.jj, attribute.recovery = PCAm.jj[, 2])
      recovery$plotPCV <- ggplot2::ggplot(data = PCV.jj, ggplot2::aes(x=item.position, y=pattern.recovery)) +
        ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(), panel.grid.minor.y = ggplot2::element_blank()) +
        ggplot2::scale_x_continuous("Until Item Position", labels = 1:MAXJ, breaks = 1:MAXJ) +
        ggplot2::scale_y_continuous("Pattern Recovery", limits = c(0,1), labels = seq(from = 0, to = 1, by = 0.10),
                                    breaks = seq(from = 0, to = 1, by = 0.10)) +
        ggplot2::geom_line() + ggplot2::geom_point(size = 2)
      recovery$plotPCAm <- ggplot2::ggplot(data = PCAm.jj, ggplot2::aes(x=item.position, y=attribute.recovery)) +
        ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(), panel.grid.minor.y = ggplot2::element_blank()) +
        ggplot2::scale_x_continuous("Until Item Position", labels = 1:MAXJ, breaks = 1:MAXJ) +
        ggplot2::scale_y_continuous("Attribute Recovery", limits = c(0,1), labels = seq(from = 0, to = 1, by = 0.10),
                                    breaks = seq(from = 0, to = 1, by = 0.10)) +
        ggplot2::geom_line() + ggplot2::geom_point(size = 2)

      J  <- nrow(cdcat.obj$specifications$Q)
      tmp <- table(unlist(lapply(cdcat.obj$est, function(x) x[[2]])))
      df.exp <- data.frame("item" = 1:J,
                           "exp" = 0)
      df.exp[as.numeric(names(tmp)), 2] <- tmp/N
      exposure <- list()
      exposure$stats <- c(summary(df.exp[, 2]), "overlap rate" = overlap.rate(MAXJ, df.exp[, 2]))
      exposure$plot <-
        ggplot2::ggplot(data=df.exp, ggplot2::aes(x=item, y=exp)) +
        ggplot2::geom_bar(stat="identity") +
        ggplot2::scale_x_continuous(name = "Item") +
        ggplot2::scale_y_continuous(name = "Exposure rate", limits = c(0, 1))
      exposure$dataplot <- df.exp

      return(list("recovery" = recovery, "item.exposure" = exposure))

    } else {

      est <- lapply(cdcat.obj$est, '[[', 1)
      N <- length(est)
      K <- ncol(cdcat.obj$specifications$Q)
      est.MAP.end <- matrix(data = as.numeric(
        unlist(lapply(lapply(lapply(est, function(x) x[nrow(x), "MAP"]), as.character), FUN = strsplit, split = ""))),
        ncol = K, nrow = N, byrow = TRUE)

      recovery <- list()
      for (kk in 1:K) {
        recovery$PCV[kk] <- mean(rowSums(est.MAP.end == alpha) >= kk)}
      recovery$PCA <- c(colMeans(est.MAP.end == alpha), mean(est.MAP.end == alpha))
      names(recovery$PCV) <- paste(1:K, sep = "/", K)
      names(recovery$PCA) <- c(paste("K", 1:K, sep = ""), "mean.PCA")

      CATlength <- list()
      data.len <- data.frame("length" = unlist(lapply(lapply(cdcat.obj$est, '[[', 2), length)),
                             "cond" = rep(1, N))
      CATlength$stats <- summary(data.len[, 1])
      CATlength$plot <- ggplot2::ggplot(data.len, ggplot2::aes(x=cond, y=length)) +
        ggplot2::geom_violin(alpha=0.4) +
        ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(), panel.grid.minor.y = ggplot2::element_blank()) +
        ggplot2::geom_dotplot(binaxis = 'y', stackdir = 'center', binwidth = 0.15) +
        ggplot2::stat_summary(fun=mean, geom="point", shape=20, size=10, color="red", fill="red") +
        ggplot2::theme(legend.position="none",
                       axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank()) +
        ggplot2::scale_fill_brewer(palette="Set3") +
        ggplot2::scale_y_continuous("CAT length", limits = c(0, (max(data.len[, 1]) + 2)),
                                    labels = seq(from = 0, to = (max(data.len[, 1]) + 2), by = 5),
                                    breaks = seq(from = 0, to = (max(data.len[, 1]) + 2), by = 5)) +
        ggplot2::scale_x_continuous("", labels = c("", "", ""), breaks = seq(0.75, 1.25, 0.25))

      J  <- nrow(cdcat.obj$specifications$Q)
      tmp <- table(unlist(lapply(cdcat.obj$est, function(x) x[[2]])))
      df.exp <- data.frame("item" = 1:J,
                           "exp" = 0)
      df.exp[as.numeric(names(tmp)), 2] <- tmp/N
      exposure <- list()
      exposure$stats <- c(summary(df.exp[, 2]), "overlap rate" = overlap.rate(CATlength$stats["Mean"], df.exp[, 2]))
      exposure$plot <-
        ggplot2::ggplot(data=df.exp, ggplot2::aes(x=item, y=exp)) +
        ggplot2::geom_bar(stat="identity") +
        ggplot2::scale_x_continuous(name = "Item") +
        ggplot2::scale_y_continuous(name = "Exposure rate", limits = c(0, 1))
      exposure$dataplot <- df.exp

      res <- list("recovery" = recovery, "item.exposure" = exposure, "CATlength" = CATlength)
      class(res) <- "cdcat.summary"

      return(res)
    }
  } else if(cdcat.obj$specifications$itemSelect == "NPS"){
    if(cdcat.obj$specifications$FIXED.LENGTH == TRUE){
      est <- lapply(cdcat.obj$est, '[[', 1)
      N <- length(est)
      K <- ncol(cdcat.obj$specifications$Q)
      est.MAP.jj <- list()
      PCV.jj <- as.data.frame(matrix(0, nrow = MAXJ, ncol = 2))
      PCAm.jj <- as.data.frame(matrix(0, nrow = MAXJ, ncol = 2))
      colnames(PCV.jj) <- c("item.position", "pattern.recovery")
      colnames(PCAm.jj) <- c("item.position", "attribute.recovery")
      for (jj in K:MAXJ) {
        est.MAP.jj[[jj]] <-
          matrix(data = as.numeric(
            unlist(lapply(lapply(lapply(est, function(x) x[jj, 4]), as.character), FUN = strsplit, split = ""))),
            ncol = K, nrow = N, byrow = TRUE)
        PCV.jj[jj, ] <- c(jj, mean(rowSums(est.MAP.jj[[jj]] == alpha) == K))
        PCAm.jj[jj, ] <- c(jj, mean(est.MAP.jj[[jj]] == alpha))
      }
      PCV.jj <- PCV.jj[-c(1:(K - 1)),]
      PCAm.jj <- PCAm.jj[-c(1:(K - 1)),]
      recovery <- list()
      recovery <- list()
      for (kk in 1:K) {
        recovery$PCV[kk] <- mean(rowSums(est.MAP.jj[[jj]] == alpha) >= kk)}
      recovery$PCA <- c(colMeans(est.MAP.jj[[jj]] == alpha), mean(est.MAP.jj[[jj]] == alpha))
      names(recovery$PCV) <- paste(1:K, sep = "/", K)
      names(recovery$PCA) <- c(paste("K", 1:K, sep = ""), "mean.PCA")

      recovery$dataplot <- cbind(PCV.jj, attribute.recovery = PCAm.jj[, 2])
      recovery$plotPCV <- ggplot2::ggplot(data = PCV.jj, ggplot2::aes(x=item.position, y=pattern.recovery)) +
        ggplot2::theme_gray() +
        ggplot2::scale_x_continuous("Until Item Position", labels = K:MAXJ, breaks = K:MAXJ) +
        ggplot2::scale_y_continuous("Pattern Recovery", limits = c(0,1), labels = seq(from = 0, to = 1, by = 0.10),
                                    breaks = seq(from = 0, to = 1, by = 0.10)) +
        ggplot2::geom_line() + ggplot2::geom_point(size = 2)
      recovery$plotPCAm <- ggplot2::ggplot(data = PCAm.jj, ggplot2::aes(x=item.position, y=attribute.recovery)) +
        ggplot2::theme_gray() +
        ggplot2::scale_x_continuous("Until Item Position", labels = K:MAXJ, breaks = K:MAXJ) +
        ggplot2::scale_y_continuous("Attribute Recovery", limits = c(0,1), labels = seq(from = 0, to = 1, by = 0.10),
                                    breaks = seq(from = 0, to = 1, by = 0.10)) +
        ggplot2::geom_line() + ggplot2::geom_point(size = 2)

      J  <- nrow(cdcat.obj$specifications$Q)
      tmp <- table(unlist(lapply(cdcat.obj$est, function(x) x[[2]])))
      df.exp <- data.frame("item" = 1:J,
                           "exp" = 0)
      df.exp[as.numeric(names(tmp)), 2] <- tmp/N
      exposure <- list()
      exposure$stats <- c(summary(df.exp[, 2]), "overlap rate" = overlap.rate(MAXJ, df.exp[, 2]))
      exposure$plot <-
        ggplot2::ggplot(data=df.exp, ggplot2::aes(x=item, y=exp)) +
        ggplot2::geom_bar(stat="identity") +
        ggplot2::scale_x_continuous(name = "Item") +
        ggplot2::scale_y_continuous(name = "Exposure rate", limits = c(0, 1))
      exposure$dataplot <- df.exp

      J  <- nrow(cdcat.obj$specifications$Q)
      tmp <- table(unlist(lapply(cdcat.obj$est, function(x) x[[2]])))
      df.exp <- data.frame("item" = 1:J,
                           "exp" = 0)
      df.exp[as.numeric(names(tmp)), 2] <- tmp/N
      exposure <- list()
      exposure$stats <- c(summary(df.exp[, 2]), "overlap rate" = overlap.rate(MAXJ, df.exp[, 2]))
      exposure$plot <-
        ggplot2::ggplot(data=df.exp, ggplot2::aes(x=item, y=exp)) +
        ggplot2::geom_bar(stat="identity") +
        ggplot2::scale_x_continuous(name = "Item") +
        ggplot2::scale_y_continuous(name = "Exposure rate", limits = c(0, 1))
      exposure$dataplot <- df.exp

      return(list("recovery" = recovery, "item.exposure" = exposure))
    } else {
      N <- length(cdcat.obj$est)
      CATlength <- list()
      data.len <- data.frame("length" = unlist(lapply(lapply(cdcat.obj$est, '[[', 2), length)), "cond" = rep(1, N))

      est <- lapply(cdcat.obj$est, '[[', 1)
      N <- length(est)
      K <- ncol(cdcat.obj$specifications$Q)
      est.MAP.end <- matrix(data = as.numeric(
        unlist(lapply(lapply(est, function(x) x[nrow(x), 4]), FUN = strsplit, split = ""))),
        ncol = K, nrow = N, byrow = TRUE)
      recovery <- GDINA::ClassRate(est.MAP.end, alpha)

      CATlength$stats <- summary(data.len[, 1])
      CATlength$plot <- ggplot2::ggplot(data.len, ggplot2::aes(x=cond, y=length)) +
        ggplot2::geom_violin(alpha=0.4) +
        ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(), panel.grid.minor.y = ggplot2::element_blank()) +
        ggplot2::geom_dotplot(binaxis = 'y', stackdir = 'center', binwidth = 0.15) +
        ggplot2::stat_summary(fun=mean, geom="point", shape=20, size=10, color="red", fill="red") +
        ggplot2::theme(legend.position="none",
                       axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank()) +
        ggplot2::scale_fill_brewer(palette="Set3") +
        ggplot2::scale_y_continuous("CAT length", limits = c(0, (max(data.len[, 1]) + 2)),
                                    labels = seq(from = 0, to = (max(data.len[, 1]) + 2), by = 5),
                                    breaks = seq(from = 0, to = (max(data.len[, 1]) + 2), by = 5)) +
        ggplot2::scale_x_continuous("", labels = c("", "", ""), breaks = seq(0.75, 1.25, 0.25))

      J  <- nrow(cdcat.obj$specifications$Q)
      tmp <- table(unlist(lapply(cdcat.obj$est, function(x) x[[2]])))
      df.exp <- data.frame("item" = 1:J,
                           "exp" = 0)
      df.exp[as.numeric(names(tmp)), 2] <- tmp/N
      exposure <- list()
      exposure$stats <- c(summary(df.exp[, 2]), "overlap rate" = overlap.rate(CATlength$stats["Mean"], df.exp[, 2]))
      exposure$plot <-
        ggplot2::ggplot(data=df.exp, ggplot2::aes(x=item, y=exp)) +
        ggplot2::geom_bar(stat="identity") +
        ggplot2::scale_x_continuous(name = "Item") +
        ggplot2::scale_y_continuous(name = "Exposure rate", limits = c(0, 1))
      exposure$dataplot <- df.exp

      res <- list("recovery" = recovery, "item.exposure" = exposure, "CATlength" = CATlength)
      class(res) <- "cdcat.summary"

      return(res)
    }
  }
}

# For gen.itembank
genQ <- function(J, K, n.id = 2, qkProp, PropWithId = TRUE, min.jk = 1, max.kcor = 0.3, seed = NULL){
  # J: number of items
  # K: number of attributes
  # qkProp: proportion of non-identity matrix items with k attributes (e.g., qkProb = c(0.1, 0.45, 0.45), 0.1 of items will have 1 attribute, 0.45 will have 2 attributes, etc.)
  # n.id: number of identity matrices inside Q-matrix
  # min.jk: minimum number of non-identity matrix items measuring each attribute
  # max.kcor: maximum permitted attribute correlation
  # seed: seed

  if(!is.null(seed)){set.seed(seed)}
  if(length(qkProp) > K){stop("length(qkProp) must be <= K")}
  if(sum(qkProp) != 1){stop("sum(qkProp) must be equal to 1")}
  pat <- GDINA::attributepattern(K)
  for(k in 1:length(qkProp)){assign(paste0("pat", k), rbind(pat[which(rowSums(pat) == k),]))}
  if((K * n.id) > J){stop("Identity matrix/matrices has/have more rows than the whole Q-matrix ((K * nI) > J)")}
  I <- matrix(rep(diag(1, K), n.id), ncol = K, byrow = TRUE)
  Q <- c()

  if(!PropWithId){
    Jq <- J - K * n.id
    Jk <- c()
    for(k in 1:length(qkProp)){Jk <- c(Jk, floor(Jq * qkProp[k]))}
    while(sum(Jk) < Jq){
      tmp <- sample(1:length(Jk), 1)
      Jk[tmp] <- Jk[tmp] + 1
    }
    for(k in 1:length(qkProp)){assign(paste0("replace", k), ifelse(Jk[k] > nrow(get(paste0("pat", k))), TRUE, FALSE))}
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
    for(k in 1:length(qkProp)){assign(paste0("replace", k), ifelse(Jk[k] > nrow(get(paste0("pat", k))), TRUE, FALSE))}
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
        k <- which(pos.k == l, arr.ind = TRUE)[,2]
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
pseudoP <- function(NPCD.obj, inv.type, w.type){
  if(!(inv.type %in% c(0, "min"))){stop("inv.type = 0, 'min'")}
  if(!(w.type %in% 1:4)){stop("w.type = 1, 2, 3, 4")}

  pattern <- NPCD.obj$pattern
  N <- nrow(NPCD.obj$alpha.est)
  K <- ncol(NPCD.obj$alpha.est)
  if(is.null(N)){
    N <- 1
    K <- length(NPCD.obj$alpha.est)
  }

  pP <- matrix(NA, nrow = N, ncol = K, dimnames = list(paste0("N", 1:N), paste0("K", 1:K)))
  for(i in 1:N){
    dist.i <- NPCD.obj$loss.matrix[,i]
    pP.ik <- c()
    for(k in 1:K){
      dist.ik <- cbind(k = pattern[,k], dist = dist.i)
      dist.ik <- dist.ik[order(dist.ik[,2]),]
      prop.dist.ik <- table.dist.ik <- table(dist.ik[,1], dist.ik[,2])[2,] / (table(dist.ik[,1], dist.ik[,2])[1,] + table(dist.ik[,1], dist.ik[,2])[2,])
      max.dist <- max(dist.ik[,2])
      if(inv.type == 0){
        min.dist <- 0
        inv.change.pre <- seq(1, -1, length.out = length(min.dist:max.dist))
        inv.change <- which(unique(dist.ik[,2]) %in% which(inv.change.pre < 0))
      } else if(inv.type == "min"){
        min.dist <- min(dist.ik[,2])
        inv.change.pre <- seq(1, -1, length.out = length(min.dist:max.dist))
        inv.change <- which(inv.change.pre < 0)[which(inv.change.pre < 0) %in% 1:length(prop.dist.ik)]
      }
      prop.dist.ik[inv.change] <- 1 - table.dist.ik[inv.change]
      opt <- length(inv.change.pre)
      switch(w.type,
             w <- 2^((opt - 1):0),
             w <- exp((opt - 1):0)
      )
      names(w) <- min.dist:max.dist
      w <- w[names(w) %in% names(prop.dist.ik)]
      pP.ik <- c(pP.ik, stats::weighted.mean(prop.dist.ik, w))
    }
    pP[i,] <- pP.ik
  }
  return(pP)
}
