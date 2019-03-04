#' Comparison of multiple cdcat objects
#'
#' This function compares differet \code{cdcat} objects in terms of classification accuracy (FIXED.LENGTH == TRUE) and/or CAT lenght FIXED.LENGTH == FALSE).
#'
#' @param cdcat.obj.l List of \code{cdcat} objects to be compared.
#' @param alpha N x K matrix with the attribute patterns to be compared to the \code{cdcat} results.
#'
#' @return \code{cdcat.comp} returns an object of class \code{cdcat.comp}.
#'
#' @export
#'
cdcat.comp <- function(cdcat.obj.l, alpha, ...)
  {
  cdmCAT.obj.sum.l <- list()
  for (cc in 1:length(cdcat.obj.l)) {
    cdmCAT.obj.sum.l[[cc]] <- cdcat.summary(cdcat.obj.l[[cc]], alpha)
  }

  MAXJ <- cdcat.obj.l[[1]]$specifications$MAXJ

  if(cdcat.obj.l[[1]]$specifications$FIXED.LENGTH == TRUE) {
  Jtot <-  MAXJ * length(cdcat.obj.l)
  data0 <- lapply(cdmCAT.obj.sum.l, function(x) x$recovery$numeric)
  datacomp <- NULL
  for (cc in 1:length(cdmCAT.obj.sum.l)) {
    model <- cdcat.obj.l[[cc]]$specifications$model
    if(cdcat.obj.l[[cc]]$specifications$GDINA.obj$control$maxitr == 0) {model <- "TRUE"}
    datacomp <- rbind(datacomp,
                      cbind(data0[[cc]], rep(model, cdcat.obj.l[[1]]$specifications$MAXJ)))
  }
  datacomp <- as.data.frame(datacomp)
  colnames(datacomp) <- c("item.position", "pattern.recovery", "attribute.recovery", "model")

  PCVplot <- ggplot(data = datacomp, aes(x=item.position, y=pattern.recovery, colour = model)) +
    theme_gray() +
    scale_x_continuous("Until Item Position", labels = 1:MAXJ, breaks = 1:MAXJ) +
    scale_y_continuous("Pattern Recovery", limits = c(0,1), labels = seq(from = 0, to = 1, by = 0.10),
                       breaks = seq(from = 0, to = 1, by = 0.10)) +
    geom_line() + geom_point()

  PCAmplot <- ggplot(data = datacomp, aes(x=item.position, y=attribute.recovery, colour = model)) +
    theme_gray() +
    scale_x_continuous("Until Item Position", labels = 1:MAXJ, breaks = 1:MAXJ) +
    scale_y_continuous("Attribute Recovery", limits = c(0,1), labels = seq(from = 0, to = 1, by = 0.10),
                       breaks = seq(from = 0, to = 1, by = 0.10)) +
    geom_line() + geom_point()

  res <- list("PCVcomp" = PCVplot, "PCAmcomp" = PCAmplot, "data" = datacomp)
  class(res) <- "cdcat.comp"

  return(res)
  }

  if(cdcat.obj.l[[1]]$specifications$FIXED.LENGTH == FALSE) {
    model <- NULL
    m <- length(cdmCAT.obj.sum.l)

    for (cc in 1:m) {
      tmp <- cdcat.obj.l[[cc]]$specifications$model
      if(cdcat.obj.l[[cc]]$specifications$GDINA.obj$control$maxitr == 0) {tmp <- "TRUE"}
      model <- c(model, tmp)
    } # end cc

    stats <- lapply(cdmCAT.obj.sum.l, function(x) x$CATlength$stats)
    plots <- lapply(cdmCAT.obj.sum.l, function(x) x$CATlength$plot)
    recovery <- lapply(cdmCAT.obj.sum.l, function(x) x$recovery)

    res.stats <- matrix(unlist(stats), nrow = m, byrow = TRUE)
    rownames(res.stats) <- model
    colnames(res.stats) <- names(stats[[1]])
    or <- order(res.stats[,"Mean"])
    res.stats <- res.stats[or, ]
    res.plots <- plot_grid(plotlist = plots[or], ncol = m, labels = model[or])
    res.recov <- matrix(unlist(recovery), nrow = m, byrow = TRUE)
    rownames(res.recov) <- model
    K <- ncol(cdcat.obj.l[[1]]$specifications$GDINA.obj$options$Q)
    colnames(res.recov) <- c("PCA", paste(1:K, "/", K, sep = ""))
    or2 <- order(res.recov[,"PCA"], decreasing = TRUE)
    res.recov <- res.recov[or2, ]

    res <- list("stats" = res.stats, "plots" = res.plots, "recovery" = res.recov)
    class(res) <- "cdcat.comp"

    return(res)
  }

}
