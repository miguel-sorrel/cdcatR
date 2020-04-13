#' Comparison of multiple cdcat objects
#'
#' This function compares different \code{cdcat} objects in terms of classification accuracy (\code{FIXED.LENGTH == TRUE}) and/or CAT lenght (\code{FIXED.LENGTH == FALSE}).
#'
#' @param cdcat.obj.l List of \code{cdcat} objects to be compared
#' @param alpha N x K matrix with the attribute patterns to be compared to the \code{cdcat} results
#' @param label labels for the \code{cdcat} objects. If NULL (by default), the models are used as labels
#'
#' @return \code{cdcat.comp} returns an object of class \code{cdcat.comp}.
#'
#' @export
#'
cdcat.comp <- function(cdcat.obj.l, alpha, label = NULL){

  if(!is.null(label)){
    if(length(label) != length(cdcat.obj.l)){stop("label and cdcat.obj.l must have the same length")}
  }

  model <- label
  cdmCAT.obj.sum.l <- list()
  for(cc in 1:length(cdcat.obj.l)){cdmCAT.obj.sum.l[[cc]] <- cdcat.summary(cdcat.obj.l[[cc]], alpha)}

  MAXJ <- cdcat.obj.l[[1]]$specifications$MAXJ

  if(cdcat.obj.l[[1]]$specifications$FIXED.LENGTH == TRUE) {
    Jtot <-  MAXJ * length(cdcat.obj.l)
    data0 <- lapply(cdmCAT.obj.sum.l, function(x) x$recovery$numeric)
    datacomp <- NULL
    for(cc in 1:length(cdmCAT.obj.sum.l)){
      if(is.null(model)){
        model.cc <- cdcat.obj.l[[cc]]$specifications$model
        if(is.null(cdcat.obj.l[[cc]]$specifications$fit$control$maxitr)){
          model.cc <- model.cc
        } else if(cdcat.obj.l[[cc]]$specifications$fit$control$maxitr == 0) {
          model.cc <- "TRUE"
        }
      } else {
        model.cc <- model[cc]
      }
      if(model.cc %in% datacomp$`rep(model.cc, nrow(data0[[cc]]))`){warning("Duplicated models. Consider providing different labels")}
      datacomp <- rbind(datacomp, cbind(data0[[cc]], rep(model.cc, nrow(data0[[cc]]))))
    }

    datacomp <- as.data.frame(datacomp)
    colnames(datacomp) <- c("item.position", "pattern.recovery", "attribute.recovery", "model")

    PCVplot <- ggplot2::ggplot(data = datacomp,
                               ggplot2::aes(x=item.position, y=pattern.recovery, colour = model)) +
      ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(),
                     panel.grid.minor.y = ggplot2::element_blank()) +
      ggplot2::scale_x_continuous("Until Item Position", labels = 1:MAXJ, breaks = 1:MAXJ) +
      ggplot2::scale_y_continuous("Pattern Recovery", limits = c(0,1), labels = seq(from = 0, to = 1, by = 0.10),
                         breaks = seq(from = 0, to = 1, by = 0.10)) +
      ggplot2::geom_line() + ggplot2::geom_point()

    PCAmplot <- ggplot2::ggplot(data = datacomp,
                                ggplot2::aes(x=item.position, y=attribute.recovery, colour = model)) +
      ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(),
                     panel.grid.minor.y = ggplot2::element_blank()) +
      ggplot2::scale_x_continuous("Until Item Position", labels = 1:MAXJ, breaks = 1:MAXJ) +
      ggplot2::scale_y_continuous("Attribute Recovery", limits = c(0,1), labels = seq(from = 0, to = 1, by = 0.10),
                         breaks = seq(from = 0, to = 1, by = 0.10)) +
      ggplot2::geom_line() + ggplot2::geom_point()

    res <- list("PCVcomp" = PCVplot, "PCAmcomp" = PCAmplot, "data" = datacomp)
    class(res) <- "cdcat.comp"

    return(res)
  }

  if(cdcat.obj.l[[1]]$specifications$FIXED.LENGTH == FALSE) {

    m <- length(cdmCAT.obj.sum.l)
    if(is.null(model)){
      for(cc in 1:m){
        tmp <- cdcat.obj.l[[cc]]$specifications$model
        if(!is.null(cdcat.obj.l[[cc]]$specifications$fit$control$maxitr)){
          if(cdcat.obj.l[[cc]]$specifications$fit$control$maxitr == 0){tmp <- "TRUE"}
        }
        model <- c(model, tmp)
      } # end cc
    }

    if(any(duplicated(model))){warning("Duplicated models. Consider providing different labels")}

    stats <- lapply(cdmCAT.obj.sum.l, function(x) x$CATlength$stats)
    plots <- lapply(cdmCAT.obj.sum.l, function(x) x$CATlength$plot)
    recovery <- lapply(cdmCAT.obj.sum.l, function(x) x$recovery)

    res.stats <- matrix(unlist(stats), nrow = m, byrow = TRUE)
    rownames(res.stats) <- model
    colnames(res.stats) <- names(stats[[1]])
    or <- order(res.stats[,"Mean"])
    res.stats <- res.stats[or, ]
    for (mm in 1:m) {
      plots[or][mm][[1]] <- plots[or][mm][[1]] + ggplot2::theme(axis.ticks.x = ggplot2::element_blank())
    }
    title <- cowplot::ggdraw() +
      cowplot::draw_label(
        paste("Models from left to right:",
              paste(model[or], sep = "", collapse = ", "), collapse = ""),
        fontface = 'bold', size = 10,
        x = 0, hjust = 0) +
      ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 7))
    plot_row <- cowplot::plot_grid(plotlist = plots[or], nrow = 1)
    res.plots <- cowplot::plot_grid(title, plot_row, ncol = 1, rel_heights = c(0.1, 1))
    res.recov <- matrix(unlist(recovery), nrow = m, byrow = TRUE)
    rownames(res.recov) <- model
    K <- ncol(cdcat.obj.l[[1]]$specifications$Q)
    colnames(res.recov) <- c("PCA", paste(1:K, "/", K, sep = ""))
    or2 <- order(res.recov[,"PCA"], decreasing = TRUE)
    res.recov <- res.recov[or2, ]

    res <- list("stats" = res.stats, "plots" = res.plots, "recovery" = res.recov)
    class(res) <- "cdcat.comp"

    return(res)
  }
}
