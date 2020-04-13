#' Plots for attribute mastery estimates
#'
#' @description This function generates a plot monitoring the attribute mastery estimates (\emph{X}: item position, \emph{Y}: mastery probability).
#' If a parametric CD-CAT has been conducted, posterior probabilites (with confident intervals) of mastering each attribute are plotted.
#' If a nonparametric CD-CAT has been conducted (and pseudo-probabilites have been computed), both nonparametric classification and pseudo-posterior probabilites (with confident intervals) of mastering each attribute are plotted.
#' Colors are used in the plots to indicate mastery (green), non-mastery (red), or uncertainty (blue).
#'
#' @param cdcat.obj An object of class \code{cdcat}
#' @param i Examinee to be plotted
#' @param k Attribute/s to be plotted. Default is NULL, which plots all attributes
#'
#' @return \code{att.plot} creates a plot.
#'
#' @export
#'
att.plot <- function(cdcat.obj, i, k = NULL){

  Q <- cdcat.obj$specifications$Q
  K <- ncol(Q)
  MAXJ <- cdcat.obj$specifications$MAXJ
  if(any(k > K)){stop("k must have values lower or equal to the number of attributes in Q-matrix")}
  if(is.null(k)){
    do.k <- 1:K
  } else {
    do.k <- k
  }

  plots <- list()
  if(cdcat.obj$specifications$itemSelect != "NPS"){
    est <- lapply(cdcat.obj$est, '[[', 1)[[i]]
    est <- as.matrix(est[, 8:(7 + K)])
    est <- apply(est, 2, as.numeric)
    est <- cbind(est, sqrt(est * (1 - est)))
    colnames(est) <- c(paste("K", 1:K, sep = ""), paste("SE.K", 1:K, sep = ""))
    est <- as.data.frame(est)

    Ji <- nrow(est)

    for(k in do.k){
      est.k <- data.frame(
        "item.position" = 1:Ji,
        "est.k" = est[, k],
        "lwr" = est[, k] - 1.96 * est[, (K + k)],
        "upr" = est[, k] + 1.96 * est[, (K + k)]
      )
      est.k$lwr[est.k$lwr <= 0] <- 0
      est.k$upr[est.k$upr >= 1] <- 1
      est.k$color <- "steelblue3"
      est.k$color[which(est.k$est.k < 0.5 & est.k$upr < 0.5)] <- "firebrick3"
      est.k$color[which(est.k$est.k > 0.5 & est.k$lwr > 0.5)] <- "seagreen3"

      plots[[k]] <- ggplot2::ggplot(data = est.k, ggplot2::aes(x = item.position, y = est.k)) +
        ggplot2::theme_gray() +
        ggplot2::scale_x_continuous("Until Item Position",
                                    labels = 1:Ji,
                                    breaks = 1:Ji) +
        ggplot2::scale_y_continuous(
          paste("K", k, sep = ""),
          limits = c(0, 1),
          labels = seq(from = 0, to = 1, by = 0.25),
          breaks = seq(from = 0, to = 1, by = 0.25)
        ) +
        ggplot2::geom_ribbon(
          data = est.k,
          ggplot2::aes(ymin = lwr, ymax = upr),
          linetype = 2,
          alpha = 0.2
        ) +
        ggplot2::geom_line() +
        ggplot2::geom_point(shape = 21, color = "black", fill = est.k$color, size = 2) +
        ggplot2::theme(legend.position = "none",
                       panel.grid.minor.x = ggplot2::element_blank(),
                       panel.grid.minor.y = ggplot2::element_blank())
    }
  } else if(cdcat.obj$specifications$itemSelect == "NPS"){
    est.pre <- lapply(cdcat.obj$est, '[[', 1)[[i]]
    est <- est.pre[K:nrow(est.pre), 4]
    est <- t(as.data.frame(lapply(strsplit(est, ""), as.numeric)))
    rownames(est) <- K:nrow(est.pre); colnames(est) <- paste0("K", 1:K)
    if(cdcat.obj$specifications$NPS.args$pseudo.prob){
      est <- cbind(est, est.pre[K:nrow(est.pre), 9:(8 + K)])
      est <- cbind(est, sqrt(est[, (K + 1):(2*K)] * (1 - est[, (K + 1):(2*K)])))
      colnames(est) <- c(paste("K", 1:K, sep = ""), paste("pP.K", 1:K, sep = ""), paste("SE.pP.K", 1:K, sep = ""))
    }

    Ji <- K - 1 + nrow(est)

    for (k in do.k) {
      est.k <- data.frame("item.position" = K:Ji, "est.k" = est[, k])
      est.k$color <- "steelblue3"
      if(cdcat.obj$specifications$NPS.args$pseudo.prob){
        est.k$pP.k <- est[, K + k]
        est.k$lwr <- est[, K + k] - 1.96 * est[, (2*K + k)]
        est.k$upr <- est[, K + k] + 1.96 * est[, (2*K + k)]
        est.k$lwr[est.k$lwr <= 0] <- 0
        est.k$upr[est.k$upr >= 1] <- 1
        est.k$color[which(est.k$est.k == 0 & est.k$upr < 0.5)] <- "firebrick3"
        est.k$color[which(est.k$est.k == 1 & est.k$lwr > 0.5)] <- "seagreen3"
      } else {
        est.k$color[which(est.k$est.k == 0)] <- "firebrick3"
        est.k$color[which(est.k$est.k == 1)] <- "seagreen3"
      }

      plots[[k]] <- ggplot2::ggplot(data = est.k, ggplot2::aes(x = item.position)) +
        ggplot2::theme_gray() +
        ggplot2::scale_x_continuous("Until Item Position",
                                    labels = K:Ji,
                                    breaks = K:Ji) +
        ggplot2::scale_y_continuous(
          paste("K", k, sep = ""),
          limits = c(0, 1),
          labels = seq(from = 0, to = 1, by = 0.5),
          breaks = seq(from = 0, to = 1, by = 0.5)
        )

      if(cdcat.obj$specifications$NPS.args$pseudo.prob){
        plots[[k]] <- plots[[k]] +
          ggplot2::geom_line(ggplot2::aes(y = pP.k), color = "gray30", linetype = "longdash") +
          ggplot2::geom_ribbon(
            data = est.k,
            ggplot2::aes(ymin = lwr, ymax = upr),
            linetype = 2,
            alpha = 0.2
          )
      }

      plots[[k]] <- plots[[k]] +
        ggplot2::geom_line(ggplot2::aes(y = est.k)) +
        ggplot2::geom_point(ggplot2::aes(y = est.k), shape = 21, color = "black", fill = est.k$color, size = 2) +
        ggplot2::theme(legend.position = "none", panel.grid.minor = ggplot2::element_blank())
    }
  }

  cowplot::plot_grid(plotlist = plots[do.k], nrow = length(do.k))
}
