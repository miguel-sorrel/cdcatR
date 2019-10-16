#' Item-level model comparison using 2LR test
#'
#' This function evaluates whether the saturated G-DINA model can be replaced by reduced CDMs without significant loss in model data fit for each item using two-step likelihood ratio test (2LR). Sorrel, de la Torre, Abad, and Olea (2017) and Ma & de la Torre (2018) can be consulted for details.
#'
#' @param GDINA.obj Calibrated item bank. An estimated model object of class GDINA
#'
#' @return \code{LR2step} returns an object of class \code{LR2step}.
#'
#' @references
#'
#' Ma, W. & de la Torre, J. (2018). Category-level model selection for the sequential G-DINA model. \emph{Journal of Educational and Behavorial Statistic, 44}, 45-77.
#'
#' Sorrel, M. A., de la Torre, J., Abad, F. J., & Olea, J. (2017). Two-step likelihood ratio test for item-level model comparison in cognitive diagnosis models. \emph{Methodology, 13}, 39-47.
#'
#' @examples
#'N <- 500
#'Q <- sim30GDINA$simQ
#'J <- nrow(Q)
#'gs <- data.frame(guess=rep(0.1,J),slip=rep(0.1,J))
#'
#'sim <- simGDINA(N,Q,gs.parm = gs,model = "DINA")
#'resGDINA <- GDINA(dat = sim$dat, Q = sim$Q,
#'                  model = "GDINA",verbose = F)
#'LR2 <- LR_2step(GDINA.obj = resGDINA);round(LR2$LR2.p, 2)
#
#' @export
#'

LR_2step <- function (GDINA.obj, ...)
{
  ## GDINA likelihood is computed from the GDINA.obj estimates
  # Notes from this function: (ui %*% theta - ci) >= 0

  Q <- GDINA.obj$options$Q
  Kjs <- rowSums(Q)
  Ks <- cumsum(2^Kjs)
  item <- which(rowSums(Q) > 1)

  J <- length(GDINA.obj$catprob.parm)
  out <- list()
  LR2 <- matrix(data = 0, nrow = J, ncol = 4)
  LR2.p <- matrix(data = 0, nrow = J, ncol = 4)
  df <- matrix(data = 0, nrow = J, ncol = 4)

  for (jj in 1:J) {

    pars.jj <-GDINA.obj$catprob.parm[[jj]]
    Rlj <- GDINA.obj$technicals$expectedCorrect
    Nlj <- GDINA.obj$technicals$expectedTotal
    Pj  <- GDINA.obj$catprob.parm[[jj]]
    Kjj <- sum(GDINA.obj$options$Q[jj, ])
    Kj <- length(Pj)
    RNjjob <-
      rbind(Rlj[jj, 1:Kj], Nlj[jj, 1:Kj], Nlj[jj, 1:Kj] - Rlj[jj, 1:Kj])
    if (ncol(RNjjob) == 4) {
      ## two-attribute items
      make.Lik.DINA.k2 <- function(vP) {
        p0 <- vP[1]
        p1 <- vP[1]
        p2 <- vP[1]
        p3 <- vP[1] + vP[2]
        l00 = (log(p0) * RNjjob[1, 1] + log(1 - p0) * (RNjjob[2, 1] - RNjjob[1, 1]))
        l10 = (log(p1) * RNjjob[1, 2] + log(1 - p1) * (RNjjob[2, 2] - RNjjob[1, 2]))
        l01 = (log(p2) * RNjjob[1, 3] + log(1 - p2) * (RNjjob[2, 3] - RNjjob[1, 3]))
        l11 = (log(p3) * RNjjob[1, 4] + log(1 - p3) * (RNjjob[2, 4] - RNjjob[1, 4]))
        l10
        l01
        l10
        l11
        L = l00 + l10 + l01 + l11
        # print(L)
        return(-L)
      }
      make.Lik.DINO.k2 <- function(vP) {
        p0 <- vP[1]
        p1 <- vP[1] + vP[2]
        p2 <- vP[1] + vP[2]
        p3 <- vP[1] + vP[2]
        l00 = (log(p0) * RNjjob[1, 1] + log(1 - p0) * (RNjjob[2, 1] - RNjjob[1, 1]))
        l10 = (log(p1) * RNjjob[1, 2] + log(1 - p1) * (RNjjob[2, 2] - RNjjob[1, 2]))
        l01 = (log(p2) * RNjjob[1, 3] + log(1 - p2) * (RNjjob[2, 3] - RNjjob[1, 3]))
        l11 = (log(p3) * RNjjob[1, 4] + log(1 - p3) * (RNjjob[2, 4] - RNjjob[1, 4]))
        l10
        l01
        l10
        l11
        L = l00 + l10 + l01 + l11
        # print(L)
        return(-L)
      }
      make.Lik.ACDM.k2 <- function(vP) {
        p0 <- vP[1]
        p1 <- vP[1] + vP[2]
        p2 <- vP[1] + vP[3]
        p3 <- vP[1] + vP[2] + vP[3]
        l00 = (log(p0) * RNjjob[1, 1] + log(1 - p0) * (RNjjob[2, 1] - RNjjob[1, 1]))
        l10 = (log(p1) * RNjjob[1, 2] + log(1 - p1) * (RNjjob[2, 2] - RNjjob[1, 2]))
        l01 = (log(p2) * RNjjob[1, 3] + log(1 - p2) * (RNjjob[2, 3] - RNjjob[1, 3]))
        l11 = (log(p3) * RNjjob[1, 4] + log(1 - p3) * (RNjjob[2, 4] - RNjjob[1, 4]))
        l10
        l01
        l10
        l11
        L = l00 + l10 + l01 + l11
        # print(L)
        return(-L)
      }
      make.Lik.GDINA.k2 <- function(pars) {
        p0 <- pars[1]
        p1 <- pars[2]
        p2 <- pars[3]
        p3 <- pars[4]
        l00 = (log(p0) * RNjjob[1, 1] + log(1 - p0) * (RNjjob[2, 1] - RNjjob[1, 1]))
        l10 = (log(p1) * RNjjob[1, 2] + log(1 - p1) * (RNjjob[2, 2] - RNjjob[1, 2]))
        l01 = (log(p2) * RNjjob[1, 3] + log(1 - p2) * (RNjjob[2, 3] - RNjjob[1, 3]))
        l11 = (log(p3) * RNjjob[1, 4] + log(1 - p3) * (RNjjob[2, 4] - RNjjob[1, 4]))
        l10
        l01
        l10
        l11
        L = l00 + l10 + l01 + l11
        # print(L)
        return(-L)
      }

      par.DINA <- constrOptim(theta = c(0.1, 0.1), f = make.Lik.DINA.k2, grad = NULL,
                              ui    = rbind(c(1, 0),
                                            c(0, 1),
                                            c(1, 1),
                                            c(-1, 0),
                                            c(0, -1),
                                            c(-1, -1)),
                              ci    = c(.0001, .0001, .0001,
                                        -.9999, -.9999 , -.9999))
      par.DINO <- constrOptim(theta = c(0.1, 0.1), f = make.Lik.DINO.k2, grad = NULL,
                              ui    = rbind(c(1, 0),
                                            c(0, 1),
                                            c(1, 1),
                                            c(-1, 0),
                                            c(0, -1),
                                            c(-1, -1)),
                              ci    = c(.0001, .0001, .0001,
                                        -.9999, -.9999 , -.9999))

      par.ACDM <- constrOptim(theta = c(0.1, 0.1, 0.1), f = make.Lik.ACDM.k2, grad = NULL,
                              ui = rbind(c(1, 0, 0),
                                         c(0, 1, 0),
                                         c(0, 0, 1),

                                         c(-1, 0, 0),
                                         c(0, -1, 0),
                                         c(0, 0, -1),

                                         c(-1, -1, -1)),

                              ci = c(.0001,.0001,.0001,
                                     -.9999,-.9999,-.9999,
                                     -.9999))


      par.GDINA.value <- make.Lik.GDINA.k2(pars = pars.jj)

    }

    if (ncol(RNjjob) == 8) {
      make.Lik.DINA.k3 <- function(vP) {

        p0 <- vP[1]
        p1 <- vP[1]
        p2 <- vP[1]
        p3 <- vP[1]
        p4 <- vP[1]
        p5 <- vP[1]
        p6 <- vP[1]
        p7 <- vP[1] + vP[2]
        l000 = (log(p0) * RNjjob[1, 1] + log(1 - p0) * (RNjjob[2, 1] - RNjjob[1, 1]))
        l100 = (log(p1) * RNjjob[1, 2] + log(1 - p1) * (RNjjob[2, 2] - RNjjob[1, 2]))
        l010 = (log(p2) * RNjjob[1, 3] + log(1 - p2) * (RNjjob[2, 3] - RNjjob[1, 3]))
        l001 = (log(p3) * RNjjob[1, 4] + log(1 - p3) * (RNjjob[2, 4] - RNjjob[1, 4]))
        l110 = (log(p4) * RNjjob[1, 5] + log(1 - p4) * (RNjjob[2, 5] - RNjjob[1, 5]))
        l101 = (log(p5) * RNjjob[1, 6] + log(1 - p5) * (RNjjob[2, 6] - RNjjob[1, 6]))
        l011 = (log(p6) * RNjjob[1, 7] + log(1 - p6) * (RNjjob[2, 7] - RNjjob[1, 7]))
        l111 = (log(p7) * RNjjob[1, 8] + log(1 - p7) * (RNjjob[2, 8] - RNjjob[1, 8]))
        L = l000 + l100 + l010 + l001 + l110 + l101 + l011 + l111
        # print(L)
        return(-L)
      }
      make.Lik.DINO.k3 <- function(vP) {

        p0 <- vP[1]
        p1 <- vP[1] + vP[2]
        p2 <- vP[1] + vP[2]
        p3 <- vP[1] + vP[2]
        p4 <- vP[1] + vP[2]
        p5 <- vP[1] + vP[2]
        p6 <- vP[1] + vP[2]
        p7 <- vP[1] + vP[2]
        l000 = (log(p0) * RNjjob[1, 1] + log(1 - p0) * (RNjjob[2, 1] - RNjjob[1, 1]))
        l100 = (log(p1) * RNjjob[1, 2] + log(1 - p1) * (RNjjob[2, 2] - RNjjob[1, 2]))
        l010 = (log(p2) * RNjjob[1, 3] + log(1 - p2) * (RNjjob[2, 3] - RNjjob[1, 3]))
        l001 = (log(p3) * RNjjob[1, 4] + log(1 - p3) * (RNjjob[2, 4] - RNjjob[1, 4]))
        l110 = (log(p4) * RNjjob[1, 5] + log(1 - p4) * (RNjjob[2, 5] - RNjjob[1, 5]))
        l101 = (log(p5) * RNjjob[1, 6] + log(1 - p5) * (RNjjob[2, 6] - RNjjob[1, 6]))
        l011 = (log(p6) * RNjjob[1, 7] + log(1 - p6) * (RNjjob[2, 7] - RNjjob[1, 7]))
        l111 = (log(p7) * RNjjob[1, 8] + log(1 - p7) * (RNjjob[2, 8] - RNjjob[1, 8]))
        L = l000 + l100 + l010 + l001 + l110 + l101 + l011 + l111

        return(-L)
      }
      make.Lik.ACDM.k3 <- function(vP) {

        p0 <- vP[1]
        p1 <- vP[1] + vP[2]
        p2 <- vP[1] + vP[3]
        p3 <- vP[1] + vP[4]
        p4 <- vP[1] + vP[2] + vP[3]
        p5 <- vP[1] + vP[2] + vP[4]
        p6 <- vP[1] + vP[3] + vP[4]
        p7 <- vP[1] + vP[2] + vP[3] + vP[4]
        l000 = (log(p0) * RNjjob[1, 1] + log(1 - p0) * (RNjjob[2, 1] - RNjjob[1, 1]))
        l100 = (log(p1) * RNjjob[1, 2] + log(1 - p1) * (RNjjob[2, 2] - RNjjob[1, 2]))
        l010 = (log(p2) * RNjjob[1, 3] + log(1 - p2) * (RNjjob[2, 3] - RNjjob[1, 3]))
        l001 = (log(p3) * RNjjob[1, 4] + log(1 - p3) * (RNjjob[2, 4] - RNjjob[1, 4]))
        l110 = (log(p4) * RNjjob[1, 5] + log(1 - p4) * (RNjjob[2, 5] - RNjjob[1, 5]))
        l101 = (log(p5) * RNjjob[1, 6] + log(1 - p5) * (RNjjob[2, 6] - RNjjob[1, 6]))
        l011 = (log(p6) * RNjjob[1, 7] + log(1 - p6) * (RNjjob[2, 7] - RNjjob[1, 7]))
        l111 = (log(p7) * RNjjob[1, 8] + log(1 - p7) * (RNjjob[2, 8] - RNjjob[1, 8]))
        L = l000 + l100 + l010 + l001 + l110 + l101 + l011 + l111

        return(-L)
      }
      make.Lik.GDINA.k3 <- function(pars) {

        p0 <- pars[1]
        p1 <- pars[2]
        p2 <- pars[3]
        p3 <- pars[4]
        p4 <- pars[5]
        p5 <- pars[6]
        p6 <- pars[7]
        p7 <- pars[8]
        l000 = (log(p0) * RNjjob[1, 1] + log(1 - p0) * (RNjjob[2, 1] - RNjjob[1, 1]))
        l100 = (log(p1) * RNjjob[1, 2] + log(1 - p1) * (RNjjob[2, 2] - RNjjob[1, 2]))
        l010 = (log(p2) * RNjjob[1, 3] + log(1 - p2) * (RNjjob[2, 3] - RNjjob[1, 3]))
        l001 = (log(p3) * RNjjob[1, 4] + log(1 - p3) * (RNjjob[2, 4] - RNjjob[1, 4]))
        l110 = (log(p4) * RNjjob[1, 5] + log(1 - p4) * (RNjjob[2, 5] - RNjjob[1, 5]))
        l101 = (log(p5) * RNjjob[1, 6] + log(1 - p5) * (RNjjob[2, 6] - RNjjob[1, 6]))
        l011 = (log(p6) * RNjjob[1, 7] + log(1 - p6) * (RNjjob[2, 7] - RNjjob[1, 7]))
        l111 = (log(p7) * RNjjob[1, 8] + log(1 - p7) * (RNjjob[2, 8] - RNjjob[1, 8]))
        L = l000 + l100 + l010 + l001 + l110 + l101 + l011 + l111

        return(-L)
      }
      par.DINA <- constrOptim(theta = c(0.1, 0.1), f = make.Lik.DINA.k3, grad = NULL,
                              ui    = rbind(c(1, 0),
                                            c(0, 1),
                                            c(1, 1),
                                            c(-1, 0),
                                            c(0, -1),
                                            c(-1, -1)),
                              ci    = c(.0001, .0001, .0001,
                                        -.9999, -.9999 , -.9999))
      par.DINO <- constrOptim(theta = c(0.1, 0.1), f = make.Lik.DINO.k3, grad = NULL,
                              ui    = rbind(c(1, 0),
                                            c(0, 1),
                                            c(1, 1),
                                            c(-1, 0),
                                            c(0, -1),
                                            c(-1, -1)),
                              ci    = c(.0001, .0001, .0001,
                                        -.9999, -.9999 , -.9999))

      par.ACDM <- constrOptim(theta = c(0.1, 0.1, 0.1, 0.1), f = make.Lik.ACDM.k3, grad = NULL,
                              ui = rbind(c(1, 0, 0, 0),
                                         c(0, 1, 0, 0),
                                         c(0, 0, 1, 0),
                                         c(0, 0, 0, 1),

                                         c(1, 0, 0, 0)*(-1),
                                         c(0, 1, 0, 0)*(-1),
                                         c(0, 0, 1, 0)*(-1),
                                         c(0, 0, 0, 1)*(-1),

                                         c(1, 1, 1, 1)*(-1)),

                              ci = c(.0001,.0001,.0001,.0001,
                                     -.9999,-.9999,-.9999,-.9999,
                                     -.9999))

      par.GDINA.value <- make.Lik.GDINA.k3(pars = pars.jj)

    }

    if (ncol(RNjjob) == 16) {
      make.Lik.DINA.k4 <- function(vP) {
        p0  <- vP[1]
        p1  <- vP[1]
        p2  <- vP[1]
        p3  <- vP[1]
        p4  <- vP[1]
        p5  <- vP[1]
        p6  <- vP[1]
        p7  <- vP[1]
        p8  <- vP[1]
        p9  <- vP[1]
        p10 <- vP[1]
        p11 <- vP[1]
        p12 <- vP[1]
        p13 <- vP[1]
        p14 <- vP[1]
        p15 <- vP[1] + vP[2]

        l0000 =  (log(p0) *  RNjjob[1, 1]  + log(1 - p0) * (RNjjob[2, 1] - RNjjob[1, 1]))
        l1000 =  (log(p1) *  RNjjob[1, 2]  + log(1 - p1) * (RNjjob[2, 2] - RNjjob[1, 2]))
        l0100 =  (log(p2) *  RNjjob[1, 3]  + log(1 - p2) * (RNjjob[2, 3] - RNjjob[1, 3]))
        l0010 =  (log(p3) *  RNjjob[1, 4]  + log(1 - p3) * (RNjjob[2, 4] - RNjjob[1, 4]))
        l0001 =  (log(p4) *  RNjjob[1, 5]  + log(1 - p4) * (RNjjob[2, 5] - RNjjob[1, 5]))
        l1100 =  (log(p5) *  RNjjob[1, 5]  + log(1 - p5) * (RNjjob[2, 6] - RNjjob[1, 6]))
        l1010  = (log(p6) *  RNjjob[1, 6]  + log(1 - p6) * (RNjjob[2, 7] - RNjjob[1, 7]))
        l1001  = (log(p7) *  RNjjob[1, 7]  + log(1 - p7) * (RNjjob[2, 8] - RNjjob[1, 8]))
        l0110  = (log(p8) *  RNjjob[1, 8]  + log(1 - p8) * (RNjjob[2, 9] - RNjjob[1, 9]))
        l0101  = (log(p9) *  RNjjob[1, 9]  + log(1 - p9) * (RNjjob[2, 10] - RNjjob[1, 10]))
        l0011  = (log(p10) * RNjjob[1, 10] + log(1 - p10) * (RNjjob[2, 11] - RNjjob[1, 11]))
        l1110  = (log(p11) * RNjjob[1, 11] + log(1 - p11) * (RNjjob[2, 12] - RNjjob[1, 12]))
        l1101  = (log(p12) * RNjjob[1, 12] + log(1 - p12) * (RNjjob[2, 13] - RNjjob[1, 13]))
        l1011  = (log(p13) * RNjjob[1, 13] + log(1 - p13) * (RNjjob[2, 14] - RNjjob[1, 14]))
        l0111  = (log(p14) * RNjjob[1, 14] + log(1 - p14) * (RNjjob[2, 15] - RNjjob[1, 15]))
        l1111  = (log(p15) * RNjjob[1, 15] + log(1 - p15) * (RNjjob[2, 16] - RNjjob[1, 16]))

        L = l0000 + l1000 + l0100 + l0010 + l0001 + l1100 + l1010 + l1001 + l0110 + l0101 + l0011 + l1110 + l1101 + l1011 + l0111 + l1111

        return(-L)
      }
      make.Lik.DINO.k4 <- function(vP) {
        p0  <- vP[1]
        p1  <- vP[1] + vP[2]
        p2  <- vP[1] + vP[2]
        p3  <- vP[1] + vP[2]
        p4  <- vP[1] + vP[2]
        p5  <- vP[1] + vP[2]
        p6  <- vP[1] + vP[2]
        p7  <- vP[1] + vP[2]
        p8  <- vP[1] + vP[2]
        p9  <- vP[1] + vP[2]
        p10 <- vP[1] + vP[2]
        p11 <- vP[1] + vP[2]
        p12 <- vP[1] + vP[2]
        p13 <- vP[1] + vP[2]
        p14 <- vP[1] + vP[2]
        p15 <- vP[1] + vP[2]

        l0000 =  (log(p0) *  RNjjob[1, 1]  + log(1 - p0) * (RNjjob[2, 1] - RNjjob[1, 1]))
        l1000 =  (log(p1) *  RNjjob[1, 2]  + log(1 - p1) * (RNjjob[2, 2] - RNjjob[1, 2]))
        l0100 =  (log(p2) *  RNjjob[1, 3]  + log(1 - p2) * (RNjjob[2, 3] - RNjjob[1, 3]))
        l0010 =  (log(p3) *  RNjjob[1, 4]  + log(1 - p3) * (RNjjob[2, 4] - RNjjob[1, 4]))
        l0001 =  (log(p4) *  RNjjob[1, 5]  + log(1 - p4) * (RNjjob[2, 5] - RNjjob[1, 5]))
        l1100 =  (log(p5) *  RNjjob[1, 5]  + log(1 - p5) * (RNjjob[2, 6] - RNjjob[1, 6]))
        l1010  = (log(p6) *  RNjjob[1, 6]  + log(1 - p6) * (RNjjob[2, 7] - RNjjob[1, 7]))
        l1001  = (log(p7) *  RNjjob[1, 7]  + log(1 - p7) * (RNjjob[2, 8] - RNjjob[1, 8]))
        l0110  = (log(p8) *  RNjjob[1, 8]  + log(1 - p8) * (RNjjob[2, 9] - RNjjob[1, 9]))
        l0101  = (log(p9) *  RNjjob[1, 9]  + log(1 - p9) * (RNjjob[2, 10] - RNjjob[1, 10]))
        l0011  = (log(p10) * RNjjob[1, 10] + log(1 - p10) * (RNjjob[2, 11] - RNjjob[1, 11]))
        l1110  = (log(p11) * RNjjob[1, 11] + log(1 - p11) * (RNjjob[2, 12] - RNjjob[1, 12]))
        l1101  = (log(p12) * RNjjob[1, 12] + log(1 - p12) * (RNjjob[2, 13] - RNjjob[1, 13]))
        l1011  = (log(p13) * RNjjob[1, 13] + log(1 - p13) * (RNjjob[2, 14] - RNjjob[1, 14]))
        l0111  = (log(p14) * RNjjob[1, 14] + log(1 - p14) * (RNjjob[2, 15] - RNjjob[1, 15]))
        l1111  = (log(p15) * RNjjob[1, 15] + log(1 - p15) * (RNjjob[2, 16] - RNjjob[1, 16]))

        L = l0000 + l1000 + l0100 + l0010 + l0001 + l1100 + l1010 + l1001 + l0110 + l0101 + l0011 + l1110 + l1101 + l1011 + l0111 + l1111
        return(-L)
      }

      make.Lik.ACDM.k4 <- function(vP) {

        p0 <-	vP[1]
        p1 <-	vP[1] +	vP[2]
        p2 <-	vP[1] +	vP[3]
        p3 <-	vP[1] +	vP[4]
        p4 <-	vP[1] +	vP[5]
        p5 <-	vP[1] +	vP[2] +	vP[3]
        p6 <-	vP[1] +	vP[2] +	vP[4]
        p7 <-	vP[1] +	vP[2] +	vP[5]
        p8 <-	vP[1] +	vP[3] +	vP[4]
        p9 <-	vP[1] +	vP[3] +	vP[5]
        p10 <-	vP[1] +	vP[4] +	vP[5]
        p11 <-	vP[1] +	vP[2] +	vP[3] +	vP[4]
        p12 <-	vP[1] +	vP[2] +	vP[3] +	vP[5]
        p13 <-	vP[1] +	vP[2] +	vP[4] +	vP[5]
        p14 <-	vP[1] +	vP[3] +	vP[4] +	vP[5]
        p15 <-	vP[1]	+ vP[2] +	vP[3] +	vP[4] +	vP[5]


        l0000 =  (log(p0) *  RNjjob[1, 1]  + log(1 - p0) * (RNjjob[2, 1] - RNjjob[1, 1]))
        l1000 =  (log(p1) *  RNjjob[1, 2]  + log(1 - p1) * (RNjjob[2, 2] - RNjjob[1, 2]))
        l0100 =  (log(p2) *  RNjjob[1, 3]  + log(1 - p2) * (RNjjob[2, 3] - RNjjob[1, 3]))
        l0010 =  (log(p3) *  RNjjob[1, 4]  + log(1 - p3) * (RNjjob[2, 4] - RNjjob[1, 4]))
        l0001 =  (log(p4) *  RNjjob[1, 5]  + log(1 - p4) * (RNjjob[2, 5] - RNjjob[1, 5]))
        l1100 =  (log(p5) *  RNjjob[1, 5]  + log(1 - p5) * (RNjjob[2, 6] - RNjjob[1, 6]))
        l1010  = (log(p6) *  RNjjob[1, 6]  + log(1 - p6) * (RNjjob[2, 7] - RNjjob[1, 7]))
        l1001  = (log(p7) *  RNjjob[1, 7]  + log(1 - p7) * (RNjjob[2, 8] - RNjjob[1, 8]))
        l0110  = (log(p8) *  RNjjob[1, 8]  + log(1 - p8) * (RNjjob[2, 9] - RNjjob[1, 9]))
        l0101  = (log(p9) *  RNjjob[1, 9]  + log(1 - p9) * (RNjjob[2, 10] - RNjjob[1, 10]))
        l0011  = (log(p10) * RNjjob[1, 10] + log(1 - p10) * (RNjjob[2, 11] - RNjjob[1, 11]))
        l1110  = (log(p11) * RNjjob[1, 11] + log(1 - p11) * (RNjjob[2, 12] - RNjjob[1, 12]))
        l1101  = (log(p12) * RNjjob[1, 12] + log(1 - p12) * (RNjjob[2, 13] - RNjjob[1, 13]))
        l1011  = (log(p13) * RNjjob[1, 13] + log(1 - p13) * (RNjjob[2, 14] - RNjjob[1, 14]))
        l0111  = (log(p14) * RNjjob[1, 14] + log(1 - p14) * (RNjjob[2, 15] - RNjjob[1, 15]))
        l1111  = (log(p15) * RNjjob[1, 15] + log(1 - p15) * (RNjjob[2, 16] - RNjjob[1, 16]))

        L = l0000 + l1000 + l0100 + l0010 + l0001 + l1100 + l1010 + l1001 + l0110 + l0101 + l0011 + l1110 + l1101 + l1011 + l0111 + l1111
        return(-L)
      }
      make.Lik.GDINA.k4 <- function(pars) {

        p0 <-	pars[1]
        p1 <-	pars[2]
        p2 <-	pars[3]
        p3 <-	pars[4]
        p4 <-	pars[5]
        p5 <-	pars[6]
        p6 <-	pars[7]
        p7 <-	pars[8]
        p8 <-	pars[9]
        p9 <-	pars[10]
        p10 <-	pars[11]
        p11 <-	pars[12]
        p12 <-	pars[13]
        p13 <-	pars[14]
        p14 <-	pars[15]
        p15 <-	pars[16]

        l0000 =  (log(p0) *  RNjjob[1, 1]  + log(1 - p0) * (RNjjob[2, 1] - RNjjob[1, 1]))
        l1000 =  (log(p1) *  RNjjob[1, 2]  + log(1 - p1) * (RNjjob[2, 2] - RNjjob[1, 2]))
        l0100 =  (log(p2) *  RNjjob[1, 3]  + log(1 - p2) * (RNjjob[2, 3] - RNjjob[1, 3]))
        l0010 =  (log(p3) *  RNjjob[1, 4]  + log(1 - p3) * (RNjjob[2, 4] - RNjjob[1, 4]))
        l0001 =  (log(p4) *  RNjjob[1, 5]  + log(1 - p4) * (RNjjob[2, 5] - RNjjob[1, 5]))
        l1100 =  (log(p5) *  RNjjob[1, 5]  + log(1 - p5) * (RNjjob[2, 6] - RNjjob[1, 6]))
        l1010  = (log(p6) *  RNjjob[1, 6]  + log(1 - p6) * (RNjjob[2, 7] - RNjjob[1, 7]))
        l1001  = (log(p7) *  RNjjob[1, 7]  + log(1 - p7) * (RNjjob[2, 8] - RNjjob[1, 8]))
        l0110  = (log(p8) *  RNjjob[1, 8]  + log(1 - p8) * (RNjjob[2, 9] - RNjjob[1, 9]))
        l0101  = (log(p9) *  RNjjob[1, 9]  + log(1 - p9) * (RNjjob[2, 10] - RNjjob[1, 10]))
        l0011  = (log(p10) * RNjjob[1, 10] + log(1 - p10) * (RNjjob[2, 11] - RNjjob[1, 11]))
        l1110  = (log(p11) * RNjjob[1, 11] + log(1 - p11) * (RNjjob[2, 12] - RNjjob[1, 12]))
        l1101  = (log(p12) * RNjjob[1, 12] + log(1 - p12) * (RNjjob[2, 13] - RNjjob[1, 13]))
        l1011  = (log(p13) * RNjjob[1, 13] + log(1 - p13) * (RNjjob[2, 14] - RNjjob[1, 14]))
        l0111  = (log(p14) * RNjjob[1, 14] + log(1 - p14) * (RNjjob[2, 15] - RNjjob[1, 15]))
        l1111  = (log(p15) * RNjjob[1, 15] + log(1 - p15) * (RNjjob[2, 16] - RNjjob[1, 16]))

        L = l0000 + l1000 + l0100 + l0010 + l0001 + l1100 + l1010 + l1001 + l0110 + l0101 + l0011 + l1110 + l1101 + l1011 + l0111 + l1111
        # print(L)
        return(-L)
      }
      par.DINA <- constrOptim(theta = c(0.1, 0.1), f = make.Lik.DINA.k4, grad = NULL,
                              ui    = rbind(c(1, 0),
                                            c(0, 1),
                                            c(1, 1),
                                            c(-1, 0),
                                            c(0, -1),
                                            c(-1, -1)),
                              ci    = c(.0001, .0001, .0001,
                                        -.9999, -.9999 , -.9999))
      par.DINO <- constrOptim(theta = c(0.1, 0.1), f = make.Lik.DINO.k4, grad = NULL,
                              ui    = rbind(c(1, 0),
                                            c(0, 1),
                                            c(1, 1),
                                            c(-1, 0),
                                            c(0, -1),
                                            c(-1, -1)),
                              ci    = c(.0001, .0001, .0001,
                                        -.9999, -.9999 , -.9999))

      par.ACDM <- constrOptim(theta = c(0.1, 0.1, 0.1, 0.1, 0.1), f = make.Lik.ACDM.k4, grad = NULL,
                              ui = rbind(c(1, 0, 0, 0, 0),
                                         c(0, 1, 0, 0, 0),
                                         c(0, 0, 1, 0, 0),
                                         c(0, 0, 0, 1, 0),
                                         c(0, 0, 0, 0, 1),

                                         c(-1, 0, 0, 0, 0),
                                         c(0, -1, 0, 0, 0),
                                         c(0, 0, -1, 0, 0),
                                         c(0, 0, 0, -1, 0),
                                         c(0, 0, 0, 0, -1),

                                         c(-1, -1, -1, -1, -1)),

                              ci = c(.0001,.0001,.0001,.0001,.0001,
                                     -.9999,-.9999,-.9999,-.9999,-.9999,
                                     -.9999))

      par.GDINA.value <- make.Lik.GDINA.k4(pars = pars.jj)

    }


    if (ncol(RNjjob) == 2) {
      LR2[jj, ]   <- cbind(jj, NA, NA, NA)
      LR2.p[jj, ] <- cbind(jj, NA, NA, NA)
      df[jj, ]    <- cbind(jj, NA, NA, NA)
    } else {

      LR2DINA <- 2 * (par.GDINA.value * (-1) - par.DINA$value * (-1))
      LR2DINO <- 2 * (par.GDINA.value * (-1) - par.DINO$value * (-1))
      LR2ACDM <- 2 * (par.GDINA.value * (-1) - par.ACDM$value * (-1))

      dfDINA  <- dfDINO <- ncol(RNjjob) - 2
      dfACDM  <- ncol(RNjjob) - (Kjj + 1)

      pDINA <-  1 - pchisq(LR2DINA, dfDINA)
      pDINO <-  1 - pchisq(LR2DINO, dfDINO)
      pACDM <-  1 - pchisq(LR2ACDM, dfACDM)

      LR2[jj, ]   <- cbind(jj, LR2DINA, LR2DINO, LR2ACDM)
      LR2.p[jj, ] <- cbind(jj, pDINA, pDINO, pACDM)
      df[jj, ]    <- cbind(jj, dfDINA, dfDINO, dfACDM)}

  }

  LR2 <- na.omit(LR2)[, -1]
  LR2.p <- na.omit(LR2.p)[, -1]
  df <- na.omit(df)[, -1]

  colnames(LR2) <- colnames(LR2.p) <- colnames(df) <- c("DINA",
                                                        "DINO", "ACDM")
  rownames(LR2) <- rownames(LR2.p) <- rownames(df) <- paste("Item",
                                                            item)
  out <- list(LR2 = LR2, LR2.p = LR2.p, df = df)


  return(out)
}

