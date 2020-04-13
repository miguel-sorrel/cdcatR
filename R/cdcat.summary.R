#' Summary information for a \code{cdcat} object
#'
#' This function provides classification accuracy (\code{FIXED.LENGTH == TRUE)} and/or CAT length (\code{FIXED.LENGTH == FALSE}) results for cdcat object.
#'
#' @param cdcat.obj An object of class \code{cdcat}
#' @param alpha N x K matrix with the attribute patterns to be compared to the \code{cdcat} results
#'
#' @return \code{cdcat.summary} returns an object of class \code{cdcat.summary}.
#'
#' @export
#'
cdcat.summary <- function(cdcat.obj, alpha){
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
        est.MAP.jj[[jj]] <- matrix(data = as.numeric(unlist(lapply(lapply(est, function(x) x[jj, "MAP"]), FUN = strsplit, split = ""))),
                                   ncol = K, nrow = N, byrow = TRUE)
        PCV.jj[jj, ] <- c(jj, GDINA::ClassRate(est.MAP.jj[[jj]], alpha)$PCV[K])
        PCAm.jj[jj, ] <- c(jj, GDINA::ClassRate(est.MAP.jj[[jj]], alpha)$PCA)
      }
      recovery <- list()
      recovery$numeric <- cbind(PCV.jj, attribute.recovery = PCAm.jj[, 2])
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

      return(list("recovery" = recovery))

    } else {
      N <- length(cdcat.obj$est)
      CATlength <- list()
      data.len <- data.frame("length" = unlist(lapply(lapply(cdcat.obj$est, '[[', 2), length)),
                             "cond" = rep(1, N))
      CATlength$stats <- summary(data.len[, 1])
      CATlength$plot <- ggplot2::ggplot(data.len, ggplot2::aes(x=cond, y=length)) +
        ggplot2::geom_violin(alpha=0.4) +
        ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(), panel.grid.minor.y = ggplot2::element_blank()) +
        ggplot2::geom_dotplot(binaxis = 'y', stackdir = 'center', binwidth = 0.15) +
        ggplot2::stat_summary(fun=mean, geom="point", shape=20, size=10, color="red", fill="red") +
        ggplot2::theme(legend.position="none") +
        ggplot2::scale_fill_brewer(palette="Set3") +
        ggplot2::scale_y_continuous("CAT length", limits = c(0, (max(data.len[, 1]) + 2)),
                                    labels = seq(from = 0, to = (max(data.len[, 1]) + 2), by = 1),
                                    breaks = seq(from = 0, to = (max(data.len[, 1]) + 2), by = 1)) +
        ggplot2::scale_x_continuous("", labels = c("", "", ""), breaks = seq(0.75, 1.25, 0.25))

      est <- lapply(cdcat.obj$est, '[[', 1)
      N <- length(est)
      K <- ncol(cdcat.obj$specifications$Q)
      est.MAP.end <- matrix(data = as.numeric(
        unlist(lapply(lapply(est, function(x) x[nrow(x), "MAP"]), FUN = strsplit, split = ""))),
        ncol = K, nrow = N, byrow = TRUE)
      recovery <- GDINA::ClassRate(est.MAP.end, alpha)

      res <- list("CATlength" = CATlength, "recovery" = recovery)
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
            unlist(lapply(lapply(est, function(x) x[jj, 4]), FUN = strsplit, split = ""))),
            ncol = K, nrow = N, byrow = TRUE)
        PCV.jj[jj, ] <- c(jj, GDINA::ClassRate(est.MAP.jj[[jj]], alpha)$PCV[K])
        PCAm.jj[jj, ] <- c(jj, GDINA::ClassRate(est.MAP.jj[[jj]], alpha)$PCA)
      }
      PCV.jj <- PCV.jj[-c(1:(K - 1)),]
      PCAm.jj <- PCAm.jj[-c(1:(K - 1)),]
      recovery <- list()
      recovery$numeric <- cbind(PCV.jj, attribute.recovery = PCAm.jj[, 2])
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

      return(list("recovery" = recovery))
    } else {
      N <- length(cdcat.obj$est)
      CATlength <- list()
      data.len <- data.frame("length" = unlist(lapply(lapply(cdcat.obj$est, '[[', 2), length)), "cond" = rep(1, N))
      CATlength$stats <- summary(data.len[, 1])
      CATlength$plot <- ggplot2::ggplot(data.len, ggplot2::aes(x=cond, y=length)) +
        ggplot2::geom_violin(alpha=0.4) +
        ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(), panel.grid.minor.y = ggplot2::element_blank()) +
        ggplot2::geom_dotplot(binaxis = 'y', stackdir = 'center', binwidth = 0.15) +
        ggplot2::stat_summary(fun=mean, geom="point", shape=20, size=10, color="red", fill="red") +
        ggplot2::theme(legend.position="none") +
        ggplot2::scale_fill_brewer(palette="Set3") +
        ggplot2::scale_y_continuous("CAT length", limits = c(0, (max(data.len[, 1]) + 2)),
                                    labels = seq(from = 0, to = (max(data.len[, 1]) + 2), by = 1),
                                    breaks = seq(from = 0, to = (max(data.len[, 1]) + 2), by = 1)) +
        ggplot2::scale_x_continuous("", labels = c("", "", ""), breaks = seq(0.75, 1.25, 0.25))

      est <- lapply(cdcat.obj$est, '[[', 1)
      N <- length(est)
      K <- ncol(cdcat.obj$specifications$Q)
      est.MAP.end <- matrix(data = as.numeric(
        unlist(lapply(lapply(est, function(x) x[nrow(x), 4]), FUN = strsplit, split = ""))),
        ncol = K, nrow = N, byrow = TRUE)
      recovery <- GDINA::ClassRate(est.MAP.end, alpha)

      res <- list("CATlength" = CATlength, "recovery" = recovery)
      class(res) <- "cdcat.summary"

      return(res)
    }
  }
}
