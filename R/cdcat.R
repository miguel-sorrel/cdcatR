#' @title Cognitively based computerized adaptive test application
#'
#' @description \code{cdcat} conducts a CD-CAT application for a given dataset. Different item selection rules can be used: the general discrimination index
#' (GDI; de la Torre & Chiu, 2016; Kaplan et al., 2015), the Jensen-Shannon divergence index (JSD; Kang et al., 2017; Minchen & de la Torre, 2016;
#' Yigit et al., 2018), the posterior-weighted Kullback-Leibler index (PWKL; Cheng, 2009), the modified PWKL index (MPWKL; Kaplan et al., 2015), the nonparametric
#' item selection method (NPS; Chang et al., 2019), the general nonparametric item selection method (GNPS; Chiu & Chang, 2021), or random selection. Fixed length or fixed precision CD-CAT can be applied. Fixed precision CD-CAT with
#' NPS and GNPS is available, by using the pseudo-posterior probability of each student mastering each attribute (experimental).
#'
#' @param fit An object of class \code{GDINA}, \code{gdina} (parametric CD-CAT), or \code{GNPC} (non-parametric CD-CAT based on GNPS). Calibrated item bank with the \code{GDINA::GDINA} (Ma & de la Torre, 2020), \code{CDM::gdina} (Robitzsch et al., 2020), or \code{cdmTools::GNPC} (Najera et al., 2022) R packages functions
#' @param dat Numeric matrix of dimensions \emph{N} number of examinees x \emph{J} number of items. Dataset to be analyzed. If \code{is.null(dat)} the data is taken data from the fit object (i.e.,  the calibration sample is used)
#' @param itemSelect Scalar character. Item selection rule: \code{GDI}, \code{JSD}, \code{MPWKL}, \code{PWKL}, \code{NPS}, \code{GNPS}, or \code{random}
#' @param MAXJ Scalar numeric. Maximum number of items to be applied regardless of the \code{FIXED.LENGTH} argument. Default is 20
#' @param FIXED.LENGTH Scalar logical. Fixed CAT-length (\code{TRUE}) or fixed-precision (\code{FALSE}) application. Default is \code{TRUE}
#' @param startRule Scalar character. Starting rule: first item is selected at random with \code{random} and first item is selected using \code{itemSelect} with \code{max}. Default is \code{random}. Seed for random is \code{NPS.args$seed} 
#' @param startK Scalar logical. Start the CAT with an identity matrix (\code{TRUE}) or not proceed with startRule from the first item (\code{FALSE}). Default is \code{FALSE}
#' @param att.prior Numeric vector of length 2^\emph{K}, where \emph{K} is the number of attributes. Prior distribution for MAP/EAP estimates. Default is uniform
#' @param initial.distr Numeric vector of length 2^\emph{K}, where \emph{K} is the number of attributes. Weighting distribution to initialize \code{itemSelect} at item position 1. Default is uniform
#' @param precision.cut Scalar numeric. Cutoff for fixed-precision (assigned pattern posterior probability > precision.cut; Hsu, Wang, & Chen, 2013). When \code{itemSelect = "NPS"} this is evaluated at the attribute level using the pseudo-posterior probabilities for each attribute (\emph{K} assigned attribute pseudo-posterior probability > precision.cut). Default is .80. A higher cutoff is recommended when \code{itemSelect = "NPS"}
#' @param NP.args A list of options when \code{itemSelect = "NPS" or "GNPS"}. \code{Q} = Q-matrix to be used in the analysis. \code{gate} = "AND" or "OR", depending on whether a conjunctive o disjunctive nonparametric CDM is used. \code{PPP} = pseudo-posterior probability of each examinee mastering each attribute (experimental). \code{w} = weight type used for computing the pseudo-posterior probability (experimental)
#' @param itemExposurecontrol  Scalar character. Item exposure control: \code{NULL} or progressive method (Barrada, Olea, Ponsoda, & Abad, 2008) with \code{"progressive"}. Default is \code{NULL}. Seed for the random component is \code{NPS.args$seed}
#' @param b Scalar numeric. Acceleration parameter for the item exposure method. Only applies if \code{itemExposurecontrol = "progressive"}. In the progressive method the first item is selected at random and the last item (i.e., \code{MAXJ}) is selected purely based on \code{itemSelect}. The rest of the items are selected combining both a random and information components. The loss of importance of the random component will be linear with \code{b = 0}, inverse exponential with \code{b < 0}, or exponential with \code{b > 0}. Thus, \code{b} allows to optimize accuracy (\code{b < 0}) or item security (\code{b > 0}). Default is 0
#' @param maxr Scalar numeric. Value should be in the range 0-1. Maximum item exposure rate that is tolerated. Default is 1. Note that for \code{maxr < 1} parallel computing cannot be implemented
#' @param itemConstraint  Scalar character. Constraints that must be satisfied by the set of items applied: \code{NULL} or attribute constraint (Henson & Douglas, 2005) with \code{"attribute"}. If \code{"attribute"} is chosen, then each attribute must be measured at least a specific number of times indicated in the \code{constraint.args$ATTRIBUTEc} argument. Default is \code{NULL}
#' @param constraint.args A list of options when \code{itemConstraint != "NULL"}. At the moment it only includes the argument \code{ATTRIBUTEc} which must be a numeric vector of length \code{ncol(Q)} indicating the minimum number of items per attribute to be administered. Default is 3 
#' @param n.cores Scalar numeric. Number of cores to be used during parallelization. Default is 2
#' @param seed Numeric vector of length 1. Some methods have a random component, so a seed is required for consistent results
#' @param print.progress Scalar logical. Prints a progress bar to the console. Default is TRUE
#'
#' @return \code{cdcat} returns an object of class \code{cdcat}.
#' \describe{
#' \item{est}{A list that contains for each examinee the mastery posterior probability estimates at each step of the CAT (\code{est.cat}) and the items applied (\code{item.usage})}
#' \item{specifications}{A list that contains all the specifications}
#' }

#' @references
#'
#' Barrada, J. R., Olea, J., Ponsoda, V., & Abad, F. J. (2008). Incorporating randomness in the Fisher information for improving item-exposure control in CATs.\emph{British Journal of Mathematical and Statistical Psychology, 61}, 493-513.
#'
#' Chang, Y.-P., Chiu, C.-Y., & Tsai, R.-C. (2019). Nonparametric CAT for CD in educational settings with small samples. \emph{Applied Psychological Measurement, 43}, 543-561.
#'
#' Cheng, Y. (2009). When cognitive diagnosis meets computerized adaptive testing: CD-CAT. \emph{Psychometrika, 74}, 619-632.
#'
#' Chiu, C. Y., & Chang, Y. P. (2021). Advances in CD-CAT: The general nonparametric item selection method. \emph{Psychometrika, 86}, 1039-1057.
#'
#' de la Torre, J., & Chiu, C. Y. (2016). General method of empirical Q-matrix validation. \emph{Psychometrika, 81}, 253-273.
#'
#' George, A. C., Robitzsch, A., Kiefer, T., Gross, J., & Uenlue, A. (2016). The R Package CDM for cognitive diagnosis models. \emph{Journal of Statistical Software, 74}, 1-24. doi:10.18637/jss.v074.i02
#'
#' Henson, R., & Douglas, J. (2005). Test construction for cognitive diagnosis. \emph{Applied Psychological Measurement, 29}, 262-277.
#'
#' Hsu, C. L., Wang, W. C., & Chen, S. Y. (2013). Variable-length computerized adaptive testing based on cognitive diagnosis models. \emph{Applied Psychological Measurement, 37}, 563-582.
#'
#' Kang, H.-A., Zhang, S., & Chang, H.-H. (2017). Dual-objective item selection criteria in cognitive diagnostic computerized adaptive testing. \emph{Journal of Educational Measurement, 54}, 165-183.
#'
#' Kaplan, M., de la Torre, J., & Barrada, J. R. (2015). New item selection methods for cognitive diagnosis computerized adaptive testing. \emph{Applied Psychological Measurement, 39}, 167-188.
#'
#' Ma, W. & de la Torre, J. (2020). GDINA: The generalized DINA model framework. R package version 2.7.9. Retrived from https://CRAN.R-project.org/package=GDINA
#'
#' Minchen, N., & de la Torre, J. (2016, July). \emph{The continuous G-DINA model and the Jensen-Shannon divergence}. Paper presented at the International Meeting of the Psychometric Society, Asheville, NC, United States.
#'
#' Nájera, P., Sorrel, M. A., & Abad, F. J. (2022). cdmTools: Useful Tools for Cognitive Diagnosis Modeling. R package version 1.0.1. https://CRAN.R-project.org/package=cdmTools
#'
#' Robitzsch, A., Kiefer, T., George, A. C., & Uenlue, A. (2020). CDM: Cognitive Diagnosis Modeling. R package version 7.5-15. https://CRAN.R-project.org/package=CDM
#'
#' Yigit, H. D., Sorrel, M. A., de la Torre, J. (2018). Computerized adaptive testing for cognitively based multiple-choice data. \emph{Applied Psychological Measurement, 43}, 388-401.
#'
#' @import foreach
#' @import doSNOW
#' @import stats
#' @import utils
#' @import cdmTools
#'
#' @examples
#' \dontshow{
#'Q <- sim180GDINA$simQ
#'dat <- sim180GDINA$simdat[1:20, ]
#'att <- sim180GDINA$simalpha[1:20, ]
#'fit <- GDINA::GDINA(dat = dat, Q = Q, verbose = 0) # GDINA package
#'
#'res.FIXJ <- cdcat(fit = fit, dat = dat, FIXED.LENGTH = TRUE,
#'                  MAXJ = 20, n.cores = 2)
#'res.FIXJ$est[[1]] # estimates for the first examinee (fixed-length)
#' }
#' \donttest{
#'######################################
#'# Example 1.                         #
#'# CD-CAT simulation for a GDINA obj  #
#'######################################
#'
#'#-----------Data----------#
#'Q <- sim180GDINA$simQ
#'K <- ncol(Q)
#'dat <- sim180GDINA$simdat
#'att <- sim180GDINA$simalpha
#'
#'#----------Model estimation----------#
#'fit <- GDINA::GDINA(dat = dat, Q = Q, verbose = 0) # GDINA package
#'#fit <- CDM::gdina(data = dat, q.matrix = Q, progress = 0) # CDM package
#'
#'#---------------CD-CAT---------------#
#'res.FIXJ <- cdcat(fit = fit, dat = dat, FIXED.LENGTH = TRUE,
#'                  MAXJ = 20, n.cores = 2)
#'res.VARJ <- cdcat(fit = fit, dat = dat, FIXED.LENGTH = FALSE,
#'                  MAXJ = 20, precision.cut = .80, n.cores = 2)
#'
#'#---------------Results--------------#
#'res.FIXJ$est[[1]] # estimates for the first examinee (fixed-length)
#'res.VARJ$est[[1]] # estimates for the first examinee (fixed-precision)
#'att.plot(cdcat.obj = res.FIXJ, i = 1) # plot for the first examinee (fixed-length)
#'att.plot(cdcat.obj = res.VARJ, i = 1) # plot  for the first examinee (fixed-precision)
#'# FIXJ summary
#'res.FIXJ.sum.real <- cdcat.summary(cdcat.obj = res.FIXJ, alpha = att) # vs. real accuracy
#'res.FIXJ.sum.real$alpha.recovery$PCV.plot
#'res.FIXJ.sum.real$item.exposure$exp.plot
#'# VARJ summary
#'res.VARJ.sum.real <- cdcat.summary(cdcat.obj = res.VARJ, alpha = att)
#'res.VARJ.sum.real$alpha.recovery$PCV
#'res.VARJ.sum.real$item.exposure$stats
#'res.VARJ.sum.real$item.exposure$length.plot
#'res.VARJ.sum.real$item.exposure$exp.plot
#'# vs. maximum observable accuracy
#'att.J <- GDINA::personparm(fit, "MAP")[, -(K+1)] # GDINA package
#'# att.J <- t(sapply(strsplit(as.character(fit$pattern$map.est), ""), as.numeric)) # CDM package
#'class.J <- GDINA::ClassRate(att, att.J) # upper-limit for accuracy
#'res.FIXJ.sum.obse <- cdcat.summary(cdcat.obj = res.FIXJ, alpha = att.J)
#'res.FIXJ.sum.obse$alpha.recovery$PCV.plot + ggplot2::geom_hline(yintercept = class.J$PCV[K],
#'                                                         color = "firebrick3")
#'res.FIXJ.sum.obse$alpha.recovery$PCA.plot + ggplot2::geom_hline(yintercept = class.J$PCA,
#'                                                         color = "firebrick3")
#'
#'######################################
#'# Example 2.                         #
#'# CD-CAT simulation for multiple     #
#'# GDINA objs and comparison of       #
#'# performance on a validation sample #
#'######################################
#'
#'#----------------Data----------------#
#'Q <- sim180combination$simQ
#'K <- ncol(Q)
#'parm <- sim180combination$specifications$item.bank$simcatprob.parm
#'dat.c <- sim180combination$simdat[,,1]
#'att.c <- sim180combination$simalpha[,,1]
#'dat.v <- sim180combination$simdat[,,2]
#'att.v <- sim180combination$simalpha[,,2]
#'
#'#-----(multiple) Model estimation----#
#'fitTRUE <- GDINA::GDINA(dat = dat.c, Q = Q, catprob.parm = parm,
#'            control = list(maxitr = 0), verbose = 0)
#'
#'fitGDINA <- GDINA::GDINA(dat = dat.c, Q = Q, verbose = 0)
#'fitDINA <- GDINA::GDINA(dat = dat.c, Q = Q, model = "DINA", verbose = 0)
#'LR2step <- LR.2step(fitGDINA)
#'models <- LR2step$models.adj.pvalues
#'fitLR2 <- GDINA::GDINA(dat = dat.c, Q = Q, model = models, verbose = 0)
#'
#'#---------------CD-CAT---------------#
#'fit.l <- list(fitTRUE, fitLR2, fitGDINA, fitDINA)
#'res.FIXJ.l <- lapply(fit.l, function(x)  cdcat(dat = dat.v,fit = x,
#'                                               FIXED.LENGTH = TRUE, n.cores = 2))
#'res.VARJ.l <- lapply(fit.l, function(x)  cdcat(dat = dat.v,fit = x,
#'                                               FIXED.LENGTH = FALSE, n.cores = 2))
#'
#'#---------------Results--------------#
#'fitbest <- GDINA::GDINA(dat = dat.v, Q = Q, catprob.parm = parm,
#'           control = list(maxitr = 1), verbose = 0)
#'fitbest.acc <- GDINA::personparm(fitbest, "MAP")[, -(K+1)]
#'class.J <- GDINA::ClassRate(att.v, fitbest.acc) # upper-limit for accuracy
#'# FIXJ comparison
#'res.FIXJ.sum <- cdcat.summary(cdcat.obj = res.FIXJ.l, alpha = att.v)
#'res.FIXJ.sum$recovery$PCVcomp + ggplot2::geom_hline(yintercept = class.J$PCV[K],
#'                                                    color = "firebrick3")
#'res.FIXJ.sum$recovery$PCAmcomp + ggplot2::geom_hline(yintercept = class.J$PCA,
#'                                                    color = "firebrick3")
#'res.FIXJ.sum$item.exposure$stats
#'res.FIXJ.sum$item.exposure$plot
#'# VARJ comparison
#'res.VARJ.sum <- cdcat.summary(cdcat.obj = res.VARJ.l, alpha = att.v)
#'res.VARJ.sum$recovery
#'res.VARJ.sum$item.exposure$stats
#'res.VARJ.sum$item.exposure$plot
#'res.VARJ.sum$CATlength$stats
#'res.VARJ.sum$CATlength$plot
#'
#'######################################
#'# Example 3.                         #
#'# Nonparametric CD-CAT for           #
#'# small-scale assessment (NPS)       #
#'######################################
#'
#'#-----------Data----------#
#'Q <- sim180DINA$simQ
#'K <- ncol(Q)
#'N <- 50
#'dat <- sim180DINA$simdat[1:N,]
#'att <- sim180DINA$simalpha[1:N,]
#'
#'#--------Nonparametric CD-CAT--------#
#'res.NPS.FIXJ <- cdcat(dat = dat, itemSelect = "NPS", FIXED.LENGTH = TRUE,
#'                      MAXJ = 25, n.cores = 2,
#'                      NP.args = list(Q = Q, gate = "AND", pseudo.prob = TRUE, w.type = 2),
#'                      seed = 12345)
#'res.NPS.VARJ <- cdcat(dat = dat, itemSelect = "NPS", FIXED.LENGTH = FALSE,
#'                      MAXJ = 25, precision.cut = 0.90, n.cores = 2,
#'                      NP.args = list(Q = Q, gate = "AND", pseudo.prob = TRUE, w.type = 2),
#'                      seed = 12345)
#'
#'#---------------Results--------------#
#'res.NPS.FIXJ$est[[1]] # estimates for the first examinee (fixed-length)
#'res.NPS.VARJ$est[[1]] # estimates for the first examinee (fixed-precision)
#'att.plot(res.NPS.FIXJ, i = 1) # plot for estimates for the first examinee (fixed-length)
#'att.plot(res.NPS.VARJ, i = 1) # plot for estimates for the first examinee (fixed-precision)
#'# FIXJ summary
#'res.NPS.FIXJ.sum.real <- cdcat.summary(cdcat.obj = res.NPS.FIXJ, alpha = att) # vs. real accuracy
#'res.NPS.FIXJ.sum.real$alpha.recovery$PCV.plot
#'res.NPS.FIXJ.sum.real$item.exposure$exp.plot
#'# VARJ summary
#'res.NPS.VARJ.sum.real <- cdcat.summary(cdcat.obj = res.NPS.VARJ, alpha = att)
#'res.NPS.VARJ.sum.real$alpha.recovery$PCV.plot
#'res.NPS.VARJ.sum.real$item.exposure$stats
#'res.NPS.VARJ.sum.real$item.exposure$length.plot
#'res.NPS.VARJ.sum.real$item.exposure$exp.plot
#'# vs. maximum observable accuracy
#'fit <- NPCD::AlphaNP(Y = dat, Q = Q, gate = "AND")
#'att.J <- fit$alpha.est
#'class.J <- GDINA::ClassRate(att, att.J) # upper-limit for accuracy
#'res.NPS.FIXJ.sum.obse <- cdcat.summary(cdcat.obj = res.NPS.FIXJ, alpha = att.J)
#'res.NPS.FIXJ.sum.obse$alpha.recovery$PCV.plot + ggplot2::geom_hline(yintercept = class.J$PCV[K],
#'                                                             color = "firebrick3")
#'res.NPS.FIXJ.sum.obse$alpha.recovery$PCA.plot + ggplot2::geom_hline(yintercept = class.J$PCA,
#'                                                             color = "firebrick3")
#'                                                             
#'######################################
#'# Example 4.                         #
#'# Nonparametric CD-CAT for           #
#'# small-scale assessment (GNPS)      #
#'######################################
#'
#'#-----------Data----------#
#'Q <- sim180DINA$simQ
#'K <- ncol(Q)
#'N <- 50
#'dat <- sim180DINA$simdat[1:N,]
#'att <- sim180DINA$simalpha[1:N,]
#'
#'#----------Model calibration----------#
#'gnpc <- cdmTools::GNPC(dat = dat, Q = Q, verbose = 0)
#'
#'#--------Nonparametric CD-CAT--------#
#'res.GNPS.FIXJ <- cdcat(fit = gnpc, dat = dat, itemSelect = "GNPS", FIXED.LENGTH = TRUE,
#'                     MAXJ = 25, n.cores = 2, 
#'                     NP.args = list(Q = Q, gate = "AND", PPP = TRUE, w.type = 2),
#'                     seed = 12345)
#'res.GNPS.VARJ <- cdcat(fit = gnpc, dat = dat, itemSelect = "GNPS", FIXED.LENGTH = FALSE,
#'                     MAXJ = 25, precision.cut = 0.90, n.cores = 2, 
#'                     NP.args = list(Q = Q, gate = "AND", PPP = TRUE, w.type = 2),
#'                       seed = 12345)
#'
#'#---------------Results--------------#
#'res.GNPS.FIXJ$est[[1]] # estimates for the first examinee (fixed-length)
#'res.GNPS.VARJ$est[[1]] # estimates for the first examinee (fixed-precision)
#'att.plot(res.GNPS.FIXJ, i = 1) # plot for estimates for the first examinee (fixed-length)
#'att.plot(res.GNPS.VARJ, i = 1) # plot for estimates for the first examinee (fixed-precision)
#'# FIXJ summary
#'res.GNPS.FIXJ.sum.real <- cdcat.summary(cdcat.obj = res.GNPS.FIXJ, alpha = att) # vs. real accuracy
#'res.GNPS.FIXJ.sum.real$alpha.recovery$PCV.plot
#'res.GNPS.FIXJ.sum.real$item.exposure$exp.plot
#'# VARJ summary
#'res.GNPS.VARJ.sum.real <- cdcat.summary(cdcat.obj = res.GNPS.VARJ, alpha = att)
#'res.GNPS.VARJ.sum.real$alpha.recovery$PCV.plot
#'res.GNPS.VARJ.sum.real$item.exposure$exp.plot
#'res.GNPS.VARJ.sum.real$item.exposure$length.plot
#'}
#' @export
#'
cdcat <- function(fit = NULL, dat = NULL, itemSelect = "GDI", MAXJ = 20,
                  FIXED.LENGTH = TRUE, startRule = "random", startK = FALSE,
                  att.prior = NULL, initial.distr = NULL, precision.cut = 0.8,
                  NP.args = list(Q = NULL, gate = NULL, PPP = TRUE, w = 2),
                  itemExposurecontrol = NULL, b = 0, maxr = 1,
                  itemConstraint = NULL, constraint.args = list(ATTRIBUTEc = NULL),
                  n.cores = 2, seed = NULL, print.progress = TRUE){
  Q <- NP.args$Q
  initial.distr.arg <- initial.distr
  if(!is.null(fit)){
    if(!is.null(fit$catprob.parm)){
      if(!is.null(fit$extra$call)){
        model <- ifelse(length(unique(fit$model)) == 1, unique(fit$model), "Combination")
        if(is.null(dat)){
          X <- fit$options$dat
        } else {
          X <- dat
        }
        if(is.null(Q)){Q <- fit$options$Q}
        LC.prob <- fit$LC.prob
      } else {
        model <- ifelse(length(unique(fit$rule)) == 1, unique(fit$rule), "Combination")
        if(is.null(dat)){
          X <- fit$data
        } else {
          X <- dat
        }
        if(is.null(Q)){Q <- fit$q.matrix}
        LC.prob <- matrix(NA, nrow = nrow(Q), ncol = 2^ncol(Q),
                          dimnames = list(paste("Item", 1:nrow(Q)), apply(GDINA::attributepattern(ncol(Q)), 1, paste, collapse = "")))
        for(j in 1:nrow(Q)){
          kj <- which(Q[j, ] == 1)
          prob.j.l <- fit$probitem[fit$probitem$itemno == j, ]$prob
          tmp <- factor(apply(GDINA::attributepattern(ncol(Q))[, kj, drop = F], 1, paste, collapse = ""))
          tmp <- factor(tmp, levels = apply(unique(GDINA::attributepattern(length(kj))), 1, paste, collapse = ""))
          for(l in 1:length(prob.j.l)){
            LC.prob[j, which(as.numeric(tmp) == l)] <- prob.j.l[l]
          }
        }
      }
      LC.prob[LC.prob > 1] <- 1
      LC.prob[LC.prob < 0] <- 0
    } else {
      model <- paste0("GNPS")
      if(is.null(dat)){dat <- fit$specifications$dat}
      X <- dat
      Q <- fit$specifications$Q
      eta.w <- fit$eta.w
    }
  } else {
    model <- paste0("NP_", NP.args$gate)
    Q <- NP.args$Q
    X <- dat
  }
  N <- nrow(X)
  q.matrix <- Q
  J <- nrow(Q)
  K <- ncol(Q)
  L <- 2^K
  # if(itemSelect == "GNPS"){
  #   pattern <- NPCD:::AlphaPermute(K)
  # } else {
  #   pattern <- GDINA::attributepattern(K)
  # }
  pattern <- GDINA::attributepattern(K)
  Lclass <- apply(pattern, MARGIN = 1, FUN = function(x){paste(x, collapse = "")})
  if(is.null(att.prior)){att.prior <- rep(1/L, L)}
  if(is.null(initial.distr.arg)){initial.distr <- matrix(data = rep(rep(1/L, L), N), nrow = N)}
  if(!is.null(initial.distr.arg)){initial.distr <- matrix(data = rep(initial.distr, N), nrow = N, byrow = TRUE)}
  if(!(class(fit) != "GDINA" | class(fit) != "gdina")){stop("fit must be an object of class 'GDINA' or 'gdina'")}
  if(!(itemSelect %in% c("GDI", "JSD", "MPWKL", "PWKL", "NPS", "GNPS", "random"))){stop("itemSelect must be one of the following: 'GDI', 'JSD', 'MPWKL', 'PWKL', 'NPS', 'GNPS', 'random'")}
  if((maxr > 1) | (maxr < 1)){stop("maxr should be in the range 0-1")}
  if(is.null(itemConstraint) == FALSE){if(is.null(constraint.args$ATTRIBUTEc)){constraint.args$ATTRIBUTEc <- 3}}
  if(itemSelect == "NPS"){
    if(maxr != 1){warning("maxr argument does not operate  with itemSelect == 'NPS'")}
    if(!is.null(itemExposurecontrol)){warning("itemExposurecontrol is set to 'NULL' with itemSelect == 'NPS'")}
    if(is.null(dat)){stop("dat required when itemSelect == 'NPS'")}
    if(is.null(Q)){stop("Q required when itemSelect == 'NPS'")}
    if(is.null(NP.args$PPP)){NP.args$PPP <- TRUE}
    if(is.null(NP.args$w)){NP.args$w <- 2} ##################################### Check
    if(is.null(NP.args$gate)){
      stop("gate == 'AND' or gate == 'OR' required in NP.args")
    } else {
      if(!(NP.args$gate %in% c("AND", "OR"))){stop("gate == 'AND' or gate == 'OR' required in NP.args")}
    }
    if(!(NP.args$PPP %in% c(TRUE, FALSE))){stop("PPP == TRUE or PPP == FALSE required in NP.args")}
    if(!FIXED.LENGTH & !NP.args$PPP){warning("FIXED.LENGHT = FALSE is not available with itemSelect == 'NPS' if NP.args$PPP == FALSE. NP.args$PPP = TRUE applied instead")}
  } else if(itemSelect == "GNPS"){
    if(is.null(fit)){stop("fit required when itemSelect != 'NPS'")}
    if(maxr != 1){warning("maxr argument does not operate  with itemSelect == 'GNPS'")}
    if(!is.null(itemExposurecontrol)){warning("itemExpo4surecontrol is set to 'NULL' with itemSelect == 'GNPS'")}
    if(is.null(NP.args$PPP)){NP.args$PPP <- TRUE}
    if(is.null(NP.args$w)){NP.args$w <- 2} ##################################### Check
    if(!(NP.args$PPP %in% c(TRUE, FALSE))){stop("PPP == TRUE or PPP == FALSE required in NP.args")}
    if(!FIXED.LENGTH & !NP.args$PPP){warning("FIXED.LENGHT = FALSE is not available with itemSelect == 'GNPS' if NP.args$PPP == FALSE. NP.args$PPP = TRUE applied instead")}
  } else {
    if(is.null(fit)){stop("fit required when itemSelect != 'NPS'")}
  }
  if(is.null(seed)){seed <- sample(1:1e+06, size = 1)}
  if(n.cores > parallel::detectCores()){stop("n.cores cannot be higher than the number of cores in the CPU")}
  if(length(att.prior) != 2^K){stop("att.prior should be of length 2^K")}
  if(sum(att.prior) != 1){stop("att.prior should add up to 1")}
  out <- list()
  item.usage <- NULL
  
  ##############################################################################
  ############################## PARAMETRIC CD-CAT #############################
  ##############################################################################
  
  if(!itemSelect %in% c("NPS", "GNPS")){
    if(maxr == 1){
      cl <- parallel::makeCluster(n.cores, type = "SOCK")
      doSNOW::registerDoSNOW(cl)
      if(print.progress == TRUE){
        pb <- txtProgressBar(max = N, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
      }
      if(print.progress == FALSE){opts <- NULL}
      out <- foreach(i = 1:N, .options.snow = opts, .export = c("GDI.M", "proggresive.f", "H", "JSD.DICO.M", "PWKL.M", "MPWKL.M"), .inorder = T) %dopar% {
        try({
          if(!is.null(seed)){set.seed(seed + i)}
          est.cat <- matrix(NA, nrow = 1, ncol = 10 + K,
                            dimnames = list(1, c("j", "qj", "xj", "ML", "nmodesML", "Lik", "MAP", "nmodesMAP", "Post", "EAP", paste("K", 1:K, sep = ""))))[-1, ]
          item.log <- NULL
          item.usage <- lapply(out, function(x) x[[2]])
          mlogPost_GDI <- initial.distr[i, ]

          GDI <- switch(itemSelect,
                        GDI = GDI.M(LC.prob = LC.prob, mlogPost_GDI),
                        JSD = JSD.DICO.M(LC.prob = LC.prob, mlogPost_GDI),
                        MPWKL = MPWKL.M(LC.prob = LC.prob, mlogPost_GDI),
                        PWKL = PWKL.M(LC.prob = LC.prob, mlogPost_GDI, point.est = sample(which(mlogPost_GDI == max(mlogPost_GDI)), 1)),
                        random = runif(J, 0, 1))
          if(!is.null(itemExposurecontrol)){
            if(itemExposurecontrol == "progressive"){GDI <- proggresive.f(jjselect = 1, MAXJ = MAXJ, b = b, GDI = GDI)}
          }
          names(GDI) <- 1:length(GDI)
          jjselect <- 1
          jjcatGDI <- IP <- NULL
          for(jjselect in 1:MAXJ){
            if(startK){
              if(jjselect == 1){
                startK.pool <- which(rowSums(Q) == 1)
                if(startRule == "random"){jjcatGDI[jjselect] <- as.numeric(sample(startK.pool, 1))}
                if(startRule == "max"){jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI[startK.pool])))}
              } else if(jjselect <= K){
                past.k <- which(colSums(Q[jjcatGDI, , drop = FALSE]) > 0)
                startK.pool <- which(rowSums(Q[, past.k, drop = FALSE]) == 0 & rowSums(Q[, -past.k, drop = FALSE]) == 1)
                jjcatGDI[jjselect] <- as.numeric(sample(startK.pool, 1))
              } else if(jjselect > K){
                jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI[eligible])))
              }
            } else {
              if(jjselect == 1){
                if(startRule == "random"){jjcatGDI[jjselect] <- as.numeric(sample(1:J, 1))}
                if(startRule == "max"){jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI)))}
              } else {
                jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI[eligible])))
              }
            }
            x.jj <- X[i, jjcatGDI[jjselect]]
            IP <- cbind(IP, x.jj)
            p.xi <- matrix(LC.prob[jjcatGDI, ], ncol = L)
            p.xi.aj1_GDI <- matrix(0, nrow = L, ncol = length(IP))
            for(l in 1L:L){
              for(jj in 1L:length(IP)){
                p.xi.aj1_GDI[l, jj] <- (p.xi[jj, l]^IP[jj]) * ((1 - p.xi[jj, l])^(1 - IP[jj]))
              }
            }
            p.xi.aj_GDI <- (t(as.matrix(apply(p.xi.aj1_GDI, 1, prod))))
            mPost_GDI <- p.xi.aj_GDI * att.prior
            msumPost_GDI <- sum(mPost_GDI)
            mlogPost_GDI <- mPost_GDI/msumPost_GDI
            p.aj_GDI <- mlogPost_GDI %*% pattern
            if(!is.null(itemConstraint)){
              Q.tmp <- as.matrix(q.matrix[jjcatGDI, ])
              if(nrow(Q.tmp) == 1){Q.tmp <- matrix(Q.tmp, nrow = 1)}
              Q.tmp.AC <- colSums(Q.tmp) < constraint.args$ATTRIBUTEc
              if(sum(Q.tmp.AC) == 1){
                cond2 <- (Q[, which(Q.tmp.AC == TRUE)] > 0)
              } else {
                cond2 <- (rowSums(Q[, which(Q.tmp.AC ==  TRUE)]) > 0)
              }
              if(sum((!(1:nrow(Q) %in% jjcatGDI)) & cond2) == 0){
                eligible <- which((!(1:nrow(Q) %in% jjcatGDI)))
              } else {
                eligible <- (1:nrow(Q))[(!(1:nrow(Q) %in% jjcatGDI)) & cond2]
              }
            } else {
              eligible <- which((!(1:nrow(Q) %in% jjcatGDI)))
            }
            
            GDI <- switch(itemSelect,
                          GDI = GDI.M(LC.prob = LC.prob, mlogPost_GDI),
                          JSD = JSD.DICO.M(LC.prob = LC.prob, mlogPost_GDI),
                          MPWKL = MPWKL.M(LC.prob = LC.prob, mlogPost_GDI),
                          PWKL = PWKL.M(LC.prob = LC.prob, mlogPost_GDI, point.est = sample(which(mlogPost_GDI == max(mlogPost_GDI)), 1)),
                          random = runif(J, 0, 1))
            if(!is.null(itemExposurecontrol)){
              if(itemExposurecontrol == "progressive"){GDI <- proggresive.f(jjselect = jjselect + 1, MAXJ = MAXJ, b = b, GDI = GDI)}
            }
            names(GDI) <- 1:length(GDI)
            est.cat <- rbind(est.cat, c(jjcatGDI[jjselect],
                                        paste(q.matrix[jjcatGDI[jjselect], ], collapse = ""),
                                        IP[jjselect], Lclass[which.max(p.xi.aj_GDI)],
                                        length(which(p.xi.aj_GDI == max(p.xi.aj_GDI))),
                                        round(max(p.xi.aj_GDI), 5), Lclass[which.max(mlogPost_GDI)],
                                        length(which(mlogPost_GDI == max(mlogPost_GDI))),
                                        round(max(mlogPost_GDI), 5),
                                        paste(as.numeric(p.aj_GDI > 0.5), collapse = ""), round(p.aj_GDI, 5)))
            jjselect <- jjselect + 1
            if((jjselect) > MAXJ){break}
            if(!FIXED.LENGTH){
              if(max(mlogPost_GDI) >= precision.cut){break}
            }
          }
          item.log <- c(item.log, jjcatGDI)
          est.cat <- as.data.frame(est.cat)
          out[[i]] <- list(est.cat = as.data.frame(est.cat), item.usage = item.log)
        })
      }
    } else {
      for(i in 1:N){
        if(isTRUE(print.progress)){cat("\r Processing examinee", i, "of", N)}
        try({
          if(!is.null(seed)){set.seed(seed + i)}
          est.cat <- matrix(NA, nrow = 1, ncol = 10 + K,
                            dimnames = list(1, c("j", "qj", "xj", "ML", "nmodesML", "Lik", "MAP", "nmodesMAP", "Post", "EAP", paste("K", 1:K, sep = ""))))[-1, ]
          item.log <- NULL
          item.usage <- lapply(out, function(x) x[[2]])
          mlogPost_GDI <- initial.distr[i, ]
          
          GDI <- switch(itemSelect,
                        GDI = GDI.M(LC.prob = LC.prob, mlogPost_GDI),
                        JSD = JSD.DICO.M(LC.prob = LC.prob, mlogPost_GDI),
                        MPWKL = MPWKL.M(LC.prob = LC.prob, mlogPost_GDI),
                        PWKL = PWKL.M(LC.prob = LC.prob, mlogPost_GDI, point.est = sample(which(mlogPost_GDI == max(mlogPost_GDI)), 1)),
                        random = runif(J, 0, 1))
          if(!is.null(itemExposurecontrol)){
            if(itemExposurecontrol == "progressive"){GDI <- proggresive.f(jjselect = 1, MAXJ = MAXJ, b = b, GDI = GDI)}
          }
          jjselect <- 1
          jjcatGDI <- IP <- NULL
          usage <- c(unlist(item.usage), jjcatGDI)
          r_i <- as.vector(table(factor(usage, levels = 1:nrow(Q)))/N)
          tmp <- rep(1, times = length(GDI))
          tmp[which((1 - r_i/maxr) < 0)] <- 0
          GDI <- tmp * GDI
          names(GDI) <- 1:length(GDI)
          for(jjselect in 1:MAXJ){
            if(startK){
              if(jjselect == 1){
                startK.pool <- which(rowSums(Q) == 1)
                if(startRule == "random"){jjcatGDI[jjselect] <- as.numeric(sample(startK.pool, 1))}
                if(startRule == "max"){jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI[startK.pool])))}
              } else if(jjselect <= K){
                past.k <- which(colSums(Q[jjcatGDI, , drop = FALSE]) > 0)
                startK.pool <- which(rowSums(Q[, past.k, drop = FALSE]) == 0 & rowSums(Q[, -past.k, drop = FALSE]) == 1)
                jjcatGDI[jjselect] <- as.numeric(sample(startK.pool, 1))
              } else if(jjselect > K){
                jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI[eligible])))
              }
            } else {
              if(jjselect == 1){
                if(startRule == "random"){jjcatGDI[jjselect] <- as.numeric(sample(1:J, 1))}
                if(startRule == "max"){jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI)))}
              } else {
                jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI[eligible])))
              }
            }
            x.jj <- X[i, jjcatGDI[jjselect]]
            IP <- cbind(IP, x.jj)
            p.xi <- matrix(LC.prob[jjcatGDI, ], ncol = L)
            p.xi.aj1_GDI <- matrix(0, nrow = L, ncol = length(IP))
            for(l in 1L:L){
              for(jj in 1L:length(IP)){
                p.xi.aj1_GDI[l, jj] <- (p.xi[jj, l]^IP[jj]) * ((1 - p.xi[jj, l])^(1 - IP[jj]))
              }
            }
            p.xi.aj_GDI <- (t(as.matrix(apply(p.xi.aj1_GDI, 1, prod))))
            mPost_GDI <- p.xi.aj_GDI * att.prior
            msumPost_GDI <- sum(mPost_GDI)
            mlogPost_GDI <- mPost_GDI/msumPost_GDI
            p.aj_GDI <- mlogPost_GDI %*% pattern
            if(is.null(itemConstraint) == FALSE){
              Q.tmp <- as.matrix(q.matrix[jjcatGDI, ])
              if(nrow(Q.tmp) == 1){Q.tmp <- matrix(Q.tmp, nrow = 1)}
              Q.tmp.AC <- colSums(Q.tmp) < constraint.args$ATTRIBUTEc
              if(sum(Q.tmp.AC) == 1){
                cond2 <- (Q[, which(Q.tmp.AC == TRUE)] > 0)
              } else {
                cond2 <- (rowSums(Q[, which(Q.tmp.AC == TRUE)]) > 0)
              }
              if(sum((!(1:nrow(Q) %in% jjcatGDI)) & cond2) == 0){
                eligible <- which((!(1:nrow(Q) %in% jjcatGDI)))
              } else {
                eligible <- (1:nrow(Q))[(!(1:nrow(Q) %in% jjcatGDI)) & cond2]
              }
            } else {
              eligible <- which((!(1:nrow(Q) %in% jjcatGDI)))
            }
            GDI <- switch(itemSelect,
                          GDI = GDI.M(LC.prob = LC.prob, mlogPost_GDI),
                          JSD = JSD.DICO.M(LC.prob = LC.prob, mlogPost_GDI),
                          MPWKL = MPWKL.M(LC.prob = LC.prob, mlogPost_GDI),
                          PWKL = PWKL.M(LC.prob = LC.prob, mlogPost_GDI, point.est = sample(which(mlogPost_GDI == max(mlogPost_GDI)), 1)),
                          random = runif(J, 0, 1))
            if(!is.null(itemExposurecontrol)){
              if(itemExposurecontrol == "progressive"){GDI <- proggresive.f(jjselect = jjselect + 1, MAXJ = MAXJ, b = b, GDI = GDI)}
            }
            usage <- c(unlist(item.usage), jjcatGDI)
            r_i <- as.vector(table(factor(usage, levels = 1:nrow(Q)))/N)
            tmp <- rep(1, times = length(GDI))
            tmp[which((1 - r_i/maxr) < 0)] <- 0
            GDI <- tmp * GDI
            names(GDI) <- 1:length(GDI)
            est.cat <- rbind(est.cat, c(jjcatGDI[jjselect],
                                        paste(q.matrix[jjcatGDI[jjselect], ], collapse = ""),
                                        IP[jjselect], Lclass[which.max(p.xi.aj_GDI)],
                                        length(which(p.xi.aj_GDI == max(p.xi.aj_GDI))),
                                        round(max(p.xi.aj_GDI), 5), Lclass[which.max(mlogPost_GDI)],
                                        length(which(mlogPost_GDI == max(mlogPost_GDI))),
                                        round(max(mlogPost_GDI), 5),
                                        paste(as.numeric(p.aj_GDI > 0.5), collapse = ""), round(p.aj_GDI, 5)))
            jjselect <- jjselect + 1
            if((jjselect) > MAXJ){break}
            if(!FIXED.LENGTH){
              if(max(mlogPost_GDI) >= precision.cut){break}
            }
          }
          item.log <- c(item.log, jjcatGDI)
          est.cat <- as.data.frame(est.cat)
          out[[i]] <- list(est.cat = as.data.frame(est.cat), item.usage = item.log)
        })
      }
    }
  } else if(itemSelect == "NPS"){
    
    ############################################################################
    ################################ NPS METHOD ################################
    ############################################################################
    
    cl <- parallel::makeCluster(n.cores, type = "SOCK")
    doSNOW::registerDoSNOW(cl)
    if(print.progress == TRUE){
      pb <- txtProgressBar(max = N, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
    }
    if(print.progress == FALSE){opts <- NULL}
    att.pat <- matrix(NA, nrow = N, ncol = K)
    exp.items <- list()
    resp.pat <- list()
    if(NP.args$gate == "AND"){gate <- rep("AND", J)}
    if(NP.args$gate == "OR"){gate <- rep("OR", J)}
    if(!FIXED.LENGTH & !NP.args$PPP){NP.args$PPP <- TRUE}
    out <- foreach(i = 1:N, .options.snow = opts, .packages = "NPCD", .export = c("PPP", "NPC.eta"), .inorder = T) %dopar% {
      try({
        if(!is.null(seed)){set.seed(seed + i)}
        if(NP.args$PPP){
          out.i <- matrix(NA, ncol = 8 + K)[-1, ]
          colnames(out.i) <- c("j", "qj", "xj", "alpha", "HD.alpha", "alpha2", "HD.alpha2", "HD.diff", paste0("PPP.K", 1:K))
        } else {
          out.i <- matrix(NA, ncol = 8)[-1, ]
          colnames(out.i) <- c("j", "qj", "xj", "alpha", "HD.alpha", "alpha2", "HD.alpha2", "HD.diff")
        }
        x <- X[i, ]
        y <- c()
        Qy <- c()
        gatey <- c()
        R0 <- as.data.frame(Q)
        gate0 <- gate
        e1 <- c(1, rep(0, K - 1))
        hit <- which(sapply(apply(R0, 1, function(x) which(all(as.double(x) == e1))), sum) == 1)
        if(length(hit) > 1){hit <- sample(hit, 1)}
        R1 <- R0[-hit, ]
        gate1 <- gate0[-hit]
        y <- c(y, x[hit])
        x <- x[-hit]
        Qy <- rbind(Qy, R0[hit, ])
        gatey <- c(gatey, gate0[hit])
        if(NP.args$PPP){
          out.i <- rbind(out.i, c(as.double(rownames(Qy)), paste(Qy[1, ], collapse = ""), y, rep(NA, 5 + K)))
        } else {
          out.i <- rbind(out.i, c(as.double(rownames(Qy)), paste(Qy[1, ], collapse = ""), y, rep(NA, 5)))
        }
        for(k in 2:K){
          hit <- NULL
          assign(paste0("e", k), c(rep(0, k - 1), 1, rep(0, K - k)))
          while(length(hit) == 0){
            qj.prev <- matrix(get(paste0("e", k)), ncol = K, byrow = T)
            for(w in 1:(k - 1)){
              qj.prev <- rbind(qj.prev, sample(c(0, 1), 1) * get(paste0("e", w)) * y[w])
            }
            qj <- colSums(qj.prev)
            qj[qj > 0] <- 1
            hit <- as.vector(which(sapply(apply(get(paste0("R", k - 1)), 1, function(x) which(all(as.double(x) == qj))), sum) == 1))
            if(length(hit) > 1){hit <- sample(hit, 1)}
          }
          assign(paste0("R", k), get(paste0("R", k - 1))[-hit, ])
          assign(paste0("gate", k), get(paste0("gate", k - 1))[-hit])
          y <- c(y, x[hit])
          x <- x[-hit]
          Qy <- rbind(Qy, get(paste0("R", k - 1))[hit, ])
          gatey <- c(gatey, get(paste0("gate", k - 1))[hit])
          if(NP.args$PPP){
            out.i <- rbind(out.i, c(as.double(rownames(Qy))[k], paste(Qy[k, ], collapse = ""), y[k], rep(NA, 5 + K)))
          } else {
            out.i <- rbind(out.i, c(as.double(rownames(Qy))[k], paste(Qy[k, ], collapse = ""), y[k], rep(NA, 5)))
          }
        }
        alphanp <- NPCD::AlphaNP(t(as.matrix(y)), Qy, gatey[1])
        pat.dist <- cbind(alphanp$pattern, alphanp$loss.matrix)
        if(K > 1){random.order <- sample(1:K, K, F)}
        for(k in random.order){
          pat.dist <- pat.dist[order(pat.dist[, k]), ]
        }
        pat.dist <- pat.dist[order(pat.dist[, K + 1]), ]
        assign(paste0("a.", K), pat.dist[order(pat.dist[, K + 1]), ][1, 1:K])
        m <- 2
        assign(paste0("a2.", K), pat.dist[order(pat.dist[, K + 1]), ][m, 1:K])
        out.i[K, 4:8] <- c(paste(get(paste0("a.", K)), collapse = ""),
                           pat.dist[order(pat.dist[, K + 1]), ][1, K + 1],
                           paste(get(paste0("a2.", K)), collapse = ""),
                           pat.dist[order(pat.dist[, K + 1]), ][2, K + 1],
                           pat.dist[order(pat.dist[, K + 1]), ][2, K + 1] - pat.dist[order(pat.dist[, K + 1]), ][1, K + 1])
        if(NP.args$PPP){
          # PPP.k <- as.vector(pseudoP(alphanp, inv.type = "min", w.type = NPS.args$w.type))
          PPP.k <- as.vector(PPP(alphanp$loss.matrix, alphanp$pattern, NP.args$w)$PPP)
          out.i[K, 9:(8 + K)] <- round(PPP.k, 5)
        }
        t <- K + 1
        while(t <= MAXJ){
          hit <- NULL
          while(is.null(hit)){
            items <- sample(1:nrow(get(paste0("R", t - 1))))
            for(j in items){
              qj <- as.double(get(paste0("R", t - 1))[j, ])
              eta.a <- NPC.eta(qj, get(paste0("a.", t - 1)), get(paste0("gate", t - 1))[j])
              eta.a2 <- NPC.eta(qj, get(paste0("a2.", t - 1)), get(paste0("gate", t - 1))[j])
              if(eta.a != eta.a2){
                hit <- as.vector(which(sapply(apply(get(paste0("R", t - 1)), 1, function(x) which(all(as.double(x) == qj))), sum) == 1))
                if(length(hit) > 1){hit <- sample(hit, 1)}
                assign(paste0("R", t), get(paste0("R", t - 1))[-hit, ])
                assign(paste0("gate", t), get(paste0("gate", t - 1))[-hit])
                y <- c(y, x[hit])
                x <- x[-hit]
                Qy <- rbind(Qy, get(paste0("R", t - 1))[hit, ])
                gatey <- c(gatey, get(paste0("gate", t - 1))[hit])
                m <- 2
                break
              } else {
                next
              }
            }
            if(is.null(hit)){
              m <- m + 1
              alphanp <- NPCD::AlphaNP(t(as.matrix(y)), Qy, gatey[1])
              pat.dist <- cbind(alphanp$pattern, alphanp$loss.matrix)
              if(K > 1){random.order <- sample(1:K, K, F)}
              for(k in random.order){pat.dist <- pat.dist[order(pat.dist[, k]), ]}
              pat.dist <- pat.dist[order(pat.dist[, K + 1]), ]
              assign(paste0("a2.", t - 1), pat.dist[order(pat.dist[, K + 1]), ][m, 1:K])
            }
          }
          alphanp <- NPCD::AlphaNP(t(as.matrix(y)), Qy, gatey[1])
          pat.dist <- cbind(alphanp$pattern, alphanp$loss.matrix)
          if(K > 1){random.order <- sample(1:K, K, F)}
          for(k in random.order){pat.dist <- pat.dist[order(pat.dist[, k]), ]}
          pat.dist <- pat.dist[order(pat.dist[, K + 1]), ]
          assign(paste0("a.", t), pat.dist[order(pat.dist[, K + 1]), ][1, 1:K])
          m <- 2
          assign(paste0("a2.", t), pat.dist[order(pat.dist[, K + 1]), ][m, 1:K])
          if(!NP.args$PPP){
            out.i <- rbind(out.i, c(as.double(rownames(Qy))[t],
                                    paste(Qy[t, ], collapse = ""),
                                    y[t],
                                    paste(get(paste0("a.", t)), collapse = ""),
                                    pat.dist[order(pat.dist[, K + 1]), ][1, K + 1],
                                    paste(get(paste0("a2.", t)), collapse = ""),
                                    pat.dist[order(pat.dist[, K + 1]), ][2, K + 1],
                                    pat.dist[order(pat.dist[, K + 1]), ][2, K + 1] - pat.dist[order(pat.dist[, K + 1]), ][1, K + 1]))
          } else {
            if(NP.args$PPP){
              # PPP.k <- as.vector(pseudoP(alphanp, inv.type = "min", w.type = NPS.args$w.type))
              PPP.k <- as.vector(PPP(alphanp$loss.matrix, alphanp$pattern, NP.args$w)$PPP)
              out.i <- rbind(out.i, c(as.double(rownames(Qy))[t],
                                      paste(Qy[t, ], collapse = ""),
                                      y[t],
                                      paste(get(paste0("a.", t)), collapse = ""),
                                      pat.dist[order(pat.dist[, K + 1]), ][1, K + 1],
                                      paste(get(paste0("a2.", t)), collapse = ""),
                                      pat.dist[order(pat.dist[, K + 1]), ][2, K + 1],
                                      pat.dist[order(pat.dist[, K + 1]), ][2, K + 1] - pat.dist[order(pat.dist[, K + 1]), ][1, K + 1],
                                      round(PPP.k, 5)))
            }
          }
          if(NP.args$PPP){
            if(!FIXED.LENGTH){
              cut.PPP.k <- PPP.k
              cut.PPP.k[cut.PPP.k < 0.5] <- 1 - cut.PPP.k[cut.PPP.k < 0.5]
              if(all(cut.PPP.k > precision.cut)){break}
            }
          }
          t <- t + 1
        }
        att.pat[i, ] <- get(paste0("a.", t - 1))
        exp.items[[i]] <- as.double(rownames(Qy))
        out.i <- as.data.frame(out.i)
        out.i$j <- as.double(as.character(out.i$j))
        out.i$qj <- as.character(out.i$qj)
        out.i$xj <- as.double(as.character(out.i$xj))
        out.i$alpha <- as.character(out.i$alpha)
        out.i$HD.alpha <- as.double(as.character(out.i$HD.alpha))
        out.i$alpha2 <- as.character(out.i$alpha2)
        out.i$HD.alpha2 <- as.double(as.character(out.i$HD.alpha2))
        out.i$HD.diff <- as.double(as.character(out.i$HD.diff))
        if(NP.args$PPP){
          for(k in 1:K){out.i[, 8 + k] <- as.double(as.character(out.i[, 8 + k]))}
        }
        out[[i]] <- list(est.cat = out.i, item.usage = exp.items[[i]])
      })
    }
  } else if(itemSelect == "GNPS"){
    
    ############################################################################
    ################################ GNPS METHOD ###############################
    ############################################################################
    
    cl <- parallel::makeCluster(n.cores, type = "SOCK")
    doSNOW::registerDoSNOW(cl)
    if(print.progress == TRUE){
      pb <- txtProgressBar(max = N, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
    }
    if(print.progress == FALSE){opts <- NULL}
    att.pat <- matrix(NA, nrow = N, ncol = K)
    exp.items <- list()
    resp.pat <- list()
    if(!FIXED.LENGTH & !NP.args$PPP){NP.args$PPP <- TRUE}
    out <- foreach(i = 1:N, .options.snow = opts, .packages = "NPCD", .export = c("PPP"), .inorder = T) %dopar% {
      try({
        if(!is.null(seed)){set.seed(seed + i)}
        # out.i <- matrix(NA, nrow = MAXJ, ncol = 8)
        # colnames(out.i) <- c("j", "qj", "xj", "alpha", "eta.alpha", "alpha2", "eta.alpha2", "DPI")
        if(NP.args$PPP){
          out.i <- matrix(NA, ncol = 8 + K)[-1, ]
          colnames(out.i) <- c("j", "qj", "xj", "alpha", "eta.alpha", "alpha2", "eta.alpha2", "DPI", paste0("PPP.K", 1:K))
        } else {
          out.i <- matrix(NA, ncol = 8)[-1, ]
          colnames(out.i) <- c("j", "qj", "xj", "alpha", "eta.alpha", "alpha2", "eta.alpha2", "DPI")
        }
        x <- X[i, ]
        y <- c()
        Qy <- c()
        R0 <- as.data.frame(Q)
        hit <- sample(1:J, size = 1)
        t <- 1
        
        while(t <= MAXJ){
          Qy <- rbind(Qy, get(paste0("R", t - 1))[hit, ])
          assign(paste0("R", t), get(paste0("R", t - 1))[-hit,])
          y <- c(y, x[hit])
          x <- x[-hit]
          if(t == 1){
            pat.dist <- cbind(pattern, (y - eta.w[as.numeric(rownames(Qy)),])^2)
          } else {
            pat.dist <- cbind(pattern, colSums((y - eta.w[as.numeric(rownames(Qy)),])^2))
          }
          if(K > 1){random.order <- sample(1:K, K, F)}
          for(k in random.order){pat.dist <- pat.dist[order(pat.dist[, k]), ]}
          pat.dist <- pat.dist[order(pat.dist[, K + 1]), ]
          assign(paste0("a.", t), pat.dist[order(pat.dist[, K + 1]), ][1, 1:K])
          m <- 2
          assign(paste0("a2.", t), pat.dist[order(pat.dist[, K + 1]), ][m, 1:K])
          if(!NP.args$PPP){
            out.i <- rbind(out.i, c(as.double(rownames(Qy))[t], paste(Qy[t, ], collapse = ""), y[t],
                                    paste(get(paste0("a.", t)), collapse = ""),
                                    round(pat.dist[order(pat.dist[, K + 1]), ][1, K + 1], 3),
                                    paste(get(paste0("a2.", t)), collapse = ""),
                                    round(pat.dist[order(pat.dist[, K + 1]), ][2, K + 1], 3),
                                    round((pat.dist[order(pat.dist[, K + 1]), ][2, K + 1] - pat.dist[order(pat.dist[, K + 1]), ][1, K + 1])^2, 3)))
          } else {
            loss <- as.matrix(colSums((y - fit$eta.w[as.numeric(rownames(Qy)),,drop = FALSE])^2))
            PPP.k <- as.vector(PPP(loss, pattern, NP.args$w)$PPP)
            out.i <- rbind(out.i, c(as.double(rownames(Qy))[t], paste(Qy[t, ], collapse = ""), y[t],
                                    paste(get(paste0("a.", t)), collapse = ""),
                                    round(pat.dist[order(pat.dist[, K + 1]), ][1, K + 1], 3),
                                    paste(get(paste0("a2.", t)), collapse = ""),
                                    round(pat.dist[order(pat.dist[, K + 1]), ][2, K + 1], 3),
                                    round((pat.dist[order(pat.dist[, K + 1]), ][2, K + 1] - pat.dist[order(pat.dist[, K + 1]), ][1, K + 1])^2, 3),
                                    round(PPP.k, 5)))
          }
          items <- as.numeric(rownames(get(paste0("R", t))))
          if(length(items) > 0){
            DPI <- (eta.w[items, paste(get(paste0("a.", t)), collapse = "")] - eta.w[items, paste(get(paste0("a2.", t)), collapse = "")])^2
            hit <- which(DPI == max(DPI))
          } else {
            break
          }
          if(length(hit) > 1){hit <- sample(hit, 1)}
          if(NP.args$PPP){
            if(!FIXED.LENGTH){
              cut.PPP.k <- PPP.k
              cut.PPP.k[cut.PPP.k < 0.5] <- 1 - cut.PPP.k[cut.PPP.k < 0.5]
              if(all(cut.PPP.k > precision.cut)){break}
            }
          }
          t <- t + 1
        }
        if(t > MAXJ){t <- MAXJ}
        att.pat[i, ] <- get(paste0("a.", t))
        exp.items[[i]] <- as.double(rownames(Qy))
        out.i <- as.data.frame(out.i)
        out.i$j <- as.double(as.character(out.i$j))
        out.i$qj <- as.character(out.i$qj)
        out.i$xj <- as.double(as.character(out.i$xj))
        out.i$alpha <- as.character(out.i$alpha)
        out.i$eta.alpha <- as.double(as.character(out.i$eta.alpha))
        out.i$alpha2 <- as.character(out.i$alpha2)
        out.i$eta.alpha2 <- as.double(as.character(out.i$eta.alpha2))
        out.i$DPI <- as.double(as.character(out.i$DPI))
        list(est.cat = out.i, item.usage = exp.items[[i]])
      })
    }
  }
  
  if((print.progress == TRUE) & (maxr == 1)){
    close(pb)
    parallel::stopCluster(cl)
  }
  res <- list()
  res$est <- out
  res$specifications <- list(fit = fit, dat = dat, Q = Q, model = model,
                             itemSelect = itemSelect, MAXJ = MAXJ, FIXED.LENGTH = FIXED.LENGTH,
                             startK = startK, startRule = startRule, att.prior = att.prior,
                             initial.distr = initial.distr.arg, precision.cut = precision.cut,
                             NP.args = NP.args, itemExposurecontrol = itemExposurecontrol,
                             b = b, maxr = maxr, itemConstraint = itemConstraint,
                             constraint.args = constraint.args, n.cores = n.cores,
                             print.progress = TRUE)
  class(res) <- "cdcat"
  return(res)
}
