#' Simulated data (155 items, a combination of DINA, DINO, and A-CDM models)
#'
#' Simulated data, Q-matrix and item parameters for a 155-item bank with 5 attributes.
#'
#' @format A list with components:
#' \describe{
#' \item{\code{simQ}}{Artificial Q-matrix. Q-matrix structure is complex (items measure up to four attributes and only 5 of them are one-attribute items)}
#' \item{\code{simcatprob.parm}}{Artificial item parameters (probability of success for each latent group). Items 1-60 are DINA items, items 61-120 are DINO items, and items 121-180 are A-CDM items}
#' \item{\code{simdat.c}}{Calibration sample dataset. Simulated responses of 500 examinees}
#' \item{\code{simalpha.c}}{Calibration sample alpha patterns. Simulated attribute patterns of 500 examinees}
#' \item{\code{simdat.v}}{Validation sample dataset. Simulated responses of 500 examinees}
#' \item{\code{simalpha.v}}{Validation sample alpha patterns. Simulated attribute patterns of 500 examinees}
#' }
#'
#'
#'
"sim155complex"
