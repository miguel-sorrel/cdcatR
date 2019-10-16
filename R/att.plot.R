#' Create plots for attribute mastery estimates
#'
#' Create plots for attribute mastery estimates (\emph{X}: item position, \emph{Y}: mastery probability).
#'
#' @param cdcat.obj An object of class \code{cdcat}.
#' @param i examinee to be plotted.
#'
#' @return \code{att.plot} creates a plot.
#'
#' @export
#'
att.plot <- function(cdcat.obj, i, ...)
  {
  Q <- cdcat.obj$specifications$GDINA.obj$options$Q
  K <- ncol(Q)
  MAXJ <- cdcat.obj$specifications$MAXJ

  plots <- list()

  est <- lapply(cdcat.obj$est, '[[', 1)[[i]]
  est <- as.matrix(est[, 6:(5 + K)])
  est <- apply(est, 2, as.numeric)
  est <- cbind(est, sqrt(est * (1 - est)))
  colnames(est) <-
    c(paste("K", 1:K, sep = ""), paste("SE.K", 1:K, sep = ""))
  est <- as.data.frame(est)

  Ji <- nrow(est)

  for (k in 1:K) {
    est.k <- data.frame(
      "item.position" = 1:Ji,
      "est.k" = est[, k],
      "lwr" = est[, k] - 1.96 * est[, (K + k)],
      "upr" = est[, k] + 1.96 * est[, (K + k)]
    )
    est.k$lwr[est.k$lwr <= 0] <- 0
    est.k$upr[est.k$upr >= 1] <- 1

    plots[[k]] <- ggplot(data = est.k, aes(x = item.position, y = est.k)) +
      theme_gray() +
      scale_x_continuous("Until Item Position",
                         labels = 1:Ji,
                         breaks = 1:Ji) +
      scale_y_continuous(
        paste("K", k, sep = ""),
        limits = c(0, 1),
        labels = seq(from = 0, to = 1, by = 0.25),
        breaks = seq(from = 0, to = 1, by = 0.25)
      ) +
      geom_line() + geom_point() +
      geom_ribbon(
        data = est.k,
        aes(ymin = lwr, ymax = upr),
        linetype = 2,
        alpha = 0.4
      )
  }

  plot_grid(plotlist = plots, nrow = K)

}


