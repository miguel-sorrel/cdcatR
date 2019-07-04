gen.Q <- function(minJ.K, complexity, ...)
  {
  # dependencies: GDINA (attributepattern)

  # minJ.K: Vector indicating the minimum number of items measuring each attribute
  # complexity: Vector indicating the maximum number of attributes being measured by an item in each column of Q

  K <- length(minJ.K)
  J <- sum(minJ.K)
  pattern <- attributepattern(K)[-1,]
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

gen.att <- function(N, K, prior = NULL, ...)
  {
  # dependencies: GDINA (attributepattern)

  # N: Number of attribute patterns to be generated
  # K: Number of attributes

  if (is.null(prior)) {
    prior<-rep(1/2^K,2^K)
  }
  prior
  att.gr <- sample(x = length(prior), size = N, replace = TRUE,
                   prob = prior)
  att <- attributepattern(K)[att.gr, ]
  return(att)
}

GDI.M <- function(LC.prob, prior) {
  GDI.est <- c(0)
  for (jj in 1:nrow(LC.prob)) {
    MPj <- sum(LC.prob[jj,]*prior)
    GDI.est[jj] <- sum(prior*(LC.prob[jj,]-MPj)^2)    }
  return(GDI.est)
}

## Jensen-Shannon Divergence (JSD): Dichotomus
H <- function(v) {
  v <- v[v > 0]
  return(sum(-v * log2(v)))
}

JSD.DICO.M <- function(LC.prob, prior) {
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

PWKL.M <- function(LC.prob, prior, point.est) {
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

MPWKL.M <- function(LC.prob, prior) {
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
