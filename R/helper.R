gen.Q <- function(minJ.K, complexity, ...){
  # dependencies: GDINA (attributepattern)

  # minJ.K: Vector indicating the minimum number of items measuring each attribute
  # complexity: Vector indicating the maximum number of attributes being measured by an item in each column of Q

  K <- length(minJ.K)
  J <- sum(minJ.K)
  pattern <- GDINA::attributepattern(K)[-1,]
  n.attributes <- apply(pattern,1,sum)

  candidates <- list()
  Q <- NULL

  for (k in 1:K) {
    candidates[[k]] <- which(n.attributes <= complexity[k] & (pattern[, k] == TRUE))
    if(length(candidates[[k]]) == 1){
      J.K0 <- rep(candidates[[k]], times = minJ.K[k])
    } else {
      J.K0 <- sample(x = c(candidates[[k]]), size = minJ.K[k], replace = TRUE)}
    J.K <- pattern[J.K0, ]
    Q <- rbind(Q, J.K)}

  return(Q)
}
gen.att <- function(N, K, prior = NULL, ...){
  # dependencies: GDINA (attributepattern)

  # N: Number of attribute patterns to be generated
  # K: Number of attributes

  if (is.null(prior)) {
    prior<-rep(1/2^K,2^K)
  }
  prior
  att.gr <- sample(x = length(prior), size = N, replace = TRUE,
                   prob = prior)
  att <- GDINA::attributepattern(K)[att.gr, ]
  return(att)
}

# Item selection procedures
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
             w <- exp((opt - 1):0),
             if(opt %% 2 != 0){
               w <- 2^(c(floor(opt/2):1, 0, 1:floor(opt/2)))
             } else {
               w <- 2^(c((floor(opt/2) - 1):0, 0:(floor(opt/2) - 1)))
             },
             if(opt %% 2 != 0){
               w <- exp(c(floor(opt/2):1, 0, 1:floor(opt/2)))
             } else {
               w <- exp(c((floor(opt/2) - 1):0, 0:(floor(opt/2) - 1)))
             }
      )
      names(w) <- min.dist:max.dist
      w <- w[names(w) %in% names(prop.dist.ik)]
      pP.ik <- c(pP.ik, weighted.mean(prop.dist.ik, w))
    }
    pP[i,] <- pP.ik
  }
  return(pP)
}
