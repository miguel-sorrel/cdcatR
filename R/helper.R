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

