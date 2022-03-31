#' Summary information for a \code{cdcat} object
#'
#' This function provides classification accuracy, item exposure, and CAT length results for \code{cdcat} object. If a list of \code{cdcat} objects is included, these objects are compared through different tables and plots.
#'
#' @param cdcat.obj An object or list of objects of class \code{cdcat}
#' @param alpha Numeric matrix of dimensions \emph{N} x \emph{K} with the reference attribute patterns used to compute attribute classification accuracy. It is expected that it will contain the true, generating alpha pattern or those estimated with the entire item bank. It is a guideline to evaluate the \code{cdcat} results
#' @param label Character vector that contains the labels for the \code{cdcat} object(s). If \code{NULL} (by default), the models are used as labels
#'
#' @return \code{cdcat.summary} returns an object of class \code{cdcat.summary}.
#' \describe{
#' \item{recovery}{A list that contains the attribute classification accuracy results calculated at the pattern- (\emph{PCV}) and attribute-levels (\emph{PCA}). Two plots monitoring these variables are provided when \code{FIXED.LENGTH = TRUE}}
#' \item{item.exposure}{A list that contains the item exposure rates results: descriptive statistics (\code{stats}) and a plot representing the item exposure rates (\code{plot}). Note that when \code{FIXED.LENGTH = FALSE} the overlap rate is calculated based on the average CAT length}
#' \item{CATlength}{If the object or list of objects of class \code{cdcat} are fixed-precision applications (i.e., \code{FIXED.LENGTH = FALSE}), this additional list is included. It contains descriptive statistics (\code{stats}) and a plot (\code{plot}) describing the CAT length}
#' }
#'
#' @export
#'
cdcat.summary <- function(cdcat.obj, alpha, label = NULL){
<<<<<<< HEAD
 if(inherits(cdcat.obj, "list")) {
=======
 if(class(cdcat.obj) == "list") {
>>>>>>> b49429ea89dd4b54520fb3c544d0a264c2a5f147
    cdcat.obj.l <- cdcat.obj
    if(!is.null(label)){
      if(length(label) != length(cdcat.obj.l)){stop("label and cdcat.obj.l must have the same length")}
    }

    if(length(unique(unlist(lapply(cdcat.obj.l, function(x) {x$specifications$FIXED.LENGTH})))) != 1){
      stop("FIXED.LENGTH must be TRUE or FALSE for all the cdcat objects included")
    }

    model <- label
    cdmCAT.obj.sum.l <- list()
    for(cc in 1:length(cdcat.obj.l)){cdmCAT.obj.sum.l[[cc]] <- cdcat.getdata(cdcat.obj.l[[cc]], alpha)}

    MAXJ <- cdcat.obj.l[[1]]$specifications$MAXJ

    if(cdcat.obj.l[[1]]$specifications$FIXED.LENGTH == TRUE) {
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

      Jtot <-  MAXJ * length(cdcat.obj.l)
      data0 <- lapply(cdmCAT.obj.sum.l, function(x) x$recovery$dataplot)
      datacomp <- NULL
      for(cc in 1:m){
        datacomp <- rbind(datacomp, cbind(data0[[cc]], rep(model[cc], nrow(data0[[cc]]))))
      }
      datacomp <- as.data.frame(datacomp)
      colnames(datacomp) <- c("item.position", "pattern.recovery", "attribute.recovery", "model")

      PCVplot <- ggplot2::ggplot(data = datacomp,
                                 ggplot2::aes(x=item.position, y=pattern.recovery, colour = model)) +
        ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(),
                       panel.grid.minor.y = ggplot2::element_blank(),
                       legend.title = ggplot2::element_blank()) +
        ggplot2::scale_x_continuous("Until Item Position", labels = 1:MAXJ, breaks = 1:MAXJ) +
        ggplot2::scale_y_continuous("Pattern Recovery", limits = c(0,1), labels = seq(from = 0, to = 1, by = 0.10),
                                    breaks = seq(from = 0, to = 1, by = 0.10)) +
        ggplot2::geom_line() + ggplot2::geom_point()

      PCAmplot <- ggplot2::ggplot(data = datacomp,
                                  ggplot2::aes(x=item.position, y=attribute.recovery, colour = model)) +
        ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(),
                       panel.grid.minor.y = ggplot2::element_blank(),
                       legend.title = ggplot2::element_blank()) +
        ggplot2::scale_x_continuous("Until Item Position", labels = 1:MAXJ, breaks = 1:MAXJ) +
        ggplot2::scale_y_continuous("Attribute Recovery", limits = c(0,1), labels = seq(from = 0, to = 1, by = 0.10),
                                    breaks = seq(from = 0, to = 1, by = 0.10)) +
        ggplot2::geom_line() + ggplot2::geom_point()

      dataexp <- dataexpplot <- NULL
      for (cc in 1:m) {
        dataexp <- rbind(dataexp, cdmCAT.obj.sum.l[[cc]]$item.exposure$stats)
        tmp <- cdmCAT.obj.sum.l[[cc]]$item.exposure$dataplot
        dataexpplot <- rbind(dataexpplot, cbind(tmp, "model" = rep(model[cc], times = nrow(tmp))))
      rownames(dataexp)[cc] <- model[cc]
      }
      dataexpplot <- as.data.frame(dataexpplot)

      plots <- lapply(cdmCAT.obj.sum.l, function(x) x$item.exposure$plot)
      for (mm in 1:m) {
        plots[mm][[1]] <- plots[mm][[1]] + ggplot2::theme(
          axis.title.x=ggplot2::element_blank(),
          axis.text.x=ggplot2::element_blank(),
          axis.ticks.x=ggplot2::element_blank(),
          axis.title.y=ggplot2::element_blank())
      }

      title_gg <- ggplot2::ggplot() +
        ggplot2::labs(title = paste("Applications from top to bottom:",
                                    paste(model, sep = "", collapse = ", "), collapse = ""),
                      subtitle = "x-axis: Item, y-axis: Item exposure rate")
      plot_row <- cowplot::plot_grid(plotlist = plots, ncol = 1)
      exp.plots <- cowplot::plot_grid(title_gg, plot_row, ncol = 1, rel_heights = c(0.15, 1))

      recovery <- list()
      tmp <- lapply(cdmCAT.obj.sum.l, function(x) x$recovery$PCV)
      recovery$PCV <- matrix(unlist(tmp), nrow = m, byrow = TRUE)
      tmp <- lapply(cdmCAT.obj.sum.l, function(x) x$recovery$PCA)
      recovery$PCA <- matrix(unlist(tmp), nrow = m, byrow = TRUE)
      colnames(recovery$PCV) <- names(cdmCAT.obj.sum.l[[1]]$recovery$PCV)
      colnames(recovery$PCA) <- names(cdmCAT.obj.sum.l[[1]]$recovery$PCA)
      rownames(recovery$PCV) <- rownames(recovery$PCA) <- model
      recovery$PCVcomp <- PCVplot
      recovery$PCAmcomp <- PCAmplot
      recovery$dataplot <- datacomp

      exposure <- list()
      exposure$stats <- dataexp
      exposure$plot <- exp.plots
      exposure$dataplot <- dataexpplot
      res <- list("recovery" = recovery, "item.exposure" = exposure)
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

      recovery <- list()
      tmp <- lapply(cdmCAT.obj.sum.l, function(x) x$recovery$PCV)
      recovery$PCV <- matrix(unlist(tmp), nrow = m, byrow = TRUE)
      tmp <- lapply(cdmCAT.obj.sum.l, function(x) x$recovery$PCA)
      recovery$PCA <- matrix(unlist(tmp), nrow = m, byrow = TRUE)
      colnames(recovery$PCV) <- names(cdmCAT.obj.sum.l[[1]]$recovery$PCV)
      colnames(recovery$PCA) <- names(cdmCAT.obj.sum.l[[1]]$recovery$PCA)
      rownames(recovery$PCV) <- rownames(recovery$PCA) <- model

      stats <- lapply(cdmCAT.obj.sum.l, function(x) x$CATlength$stats)
      plots <- lapply(cdmCAT.obj.sum.l, function(x) x$CATlength$plot)
      res.stats <- matrix(unlist(stats), nrow = m, byrow = TRUE)
      rownames(res.stats) <- model
      colnames(res.stats) <- names(stats[[1]])
      for (mm in 1:m) {
        plots[mm][[1]] <- plots[mm][[1]] + ggplot2::theme(axis.ticks.x = ggplot2::element_blank())
      }
      title <- cowplot::ggdraw() +
        cowplot::draw_label(
          paste("Applications from left to right:",
                paste(model, sep = "", collapse = ", "), collapse = ""),
          fontface = 'bold', size = 10,
          x = 0, hjust = 0) +
        ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 7))
      plot_row <- cowplot::plot_grid(plotlist = plots, nrow = 1)
      res.plots <- cowplot::plot_grid(title, plot_row, ncol = 1, rel_heights = c(0.1, 1))

      dataexp <- dataexpplot <- NULL
      for (cc in 1:m) {
        dataexp <- rbind(dataexp, cdmCAT.obj.sum.l[[cc]]$item.exposure$stats)
        tmp <- cdmCAT.obj.sum.l[[cc]]$item.exposure$dataplot
        dataexpplot <- rbind(dataexpplot, cbind(tmp, "model" = rep(model[cc], times = nrow(tmp))))
        rownames(dataexp)[cc] <- model[cc]
      }
      dataexpplot <- as.data.frame(dataexpplot)

      plots <- lapply(cdmCAT.obj.sum.l, function(x) x$item.exposure$plot)
      for (mm in 1:m) {
        plots[mm][[1]] <- plots[mm][[1]] + ggplot2::theme(
          axis.title.x=ggplot2::element_blank(),
          axis.text.x=ggplot2::element_blank(),
          axis.ticks.x=ggplot2::element_blank(),
          axis.title.y=ggplot2::element_blank())
      }

      title_gg <- ggplot2::ggplot() +
        ggplot2::labs(title = paste("Applications from top to bottom:",
                                    paste(model, sep = "", collapse = ", "), collapse = ""),
                      subtitle = "x-axis: Item, y-axis: Item exposure rate")
      plot_row <- cowplot::plot_grid(plotlist = plots, ncol = 1)
      exp.plots <- cowplot::plot_grid(title_gg, plot_row, ncol = 1, rel_heights = c(0.15, 1))

      CATlength <- list()
      CATlength$stats <- res.stats
      CATlength$plot <- res.plots

      exposure <- list()
      exposure$stats <- dataexp
      exposure$plot <- exp.plots
      exposure$dataplot <- dataexpplot

      res <- list("recovery" = recovery, "item.exposure" = exposure, "CATlength" = CATlength)
    }
  } else {
    res <- cdcat.getdata(cdcat.obj, alpha)
  }
  class(res) <- "cdcat.summary"
  return(res)
}

