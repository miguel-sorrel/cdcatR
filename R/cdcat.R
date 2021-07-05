#' @title Cognitively based computerized adaptive test application
#'
#' @description \code{cdcat} conducts a CD-CAT application for a given dataset. Different item selection rules can be used: the general discrimination index
#' (GDI; de la Torre & Chiu, 2016; Kaplan et al., 2015), the Jensen-Shannon divergence index (JSD; Kang et al., 2017; Minchen & de la Torre, 2016;
#' Yigit et al., 2018), the posterior-weighted Kullback-Leibler index (PWKL; Cheng, 2009), the modified PWKL index (MPWKL; Kaplan et al., 2015), the nonparametric
#' item selection method (NPS; Chang et al., 2019), or random selection. Fixed length or fixed precision CD-CAT can be applied. Fixed precision CD-CAT with
#' NPS is available, by using the pseudo-posterior probability of each student mastering each attribute (experimental).
#'
#' @param fit An object of class \code{GDINA} or \code{gdina}. Calibrated item bank with the \code{GDINA::GDINA} (Ma & de la Torre, 2020) or \code{CDM::gdina} (Robitzsch et al., 2020) R packages functions
#' @param dat Numeric matrix of dimensions \emph{N} number of examinees x \emph{J} number of items. Dataset to be analyzed. If \code{is.null(dat)} the data is taken data from the fit object (i.e.,  the calibration sample is used)
#' @param itemSelect Scalar character. Item selection rule: \code{GDI}, \code{JSD}, \code{MPWKL}, \code{PWKL}, \code{NPS}, or \code{random}
#' @param MAXJ Scalar numeric. Maximum number of items to be applied regardless of the \code{FIXED.LENGTH} argument. Default is 20
#' @param FIXED.LENGTH Scalar logical. Fixed CAT-length (\code{TRUE}) or fixed-precision (\code{FALSE}) application. Default is \code{TRUE}
#' @param startRule Scalar character. Starting rule: first item is selected at random with \code{random} and first item is selected using \code{itemSelect} with \code{max}. Default is \code{random}. Seed for random is \code{NPS.args$seed}    
#' @param startK Scalar logical. Start the CAT with an identity matrix (\code{TRUE}) or not proceed with startRule from the first item (\code{FALSE}). Default is \code{FALSE}
#' @param att.prior Numeric vector of length 2^\emph{K}, where \emph{K} is the number of attributes. Prior distribution for MAP/EAP estimates. Default is uniform
#' @param initial.distr Numeric vector of length 2^\emph{K}, where \emph{K} is the number of attributes. Weighting distribution to initialize \code{itemSelect} at item position 1. Default is uniform
#' @param precision.cut Scalar numeric. Cutoff for fixed-precision (assigned pattern posterior probability > precision.cut; Hsu, Wang, & Chen, 2013). When \code{itemSelect = "NPS"} this is evaluated at the attribute level using the pseudo-posterior probabilities for each attribute (\emph{K} assigned attribute pseudo-posterior probability > precision.cut). Default is .80. A higher cutoff is recommended when \code{itemSelect = "NPS"}
#' @param NPS.args A list of options when \code{itemSelect = "NPS"}. \code{Q} = Q-matrix to be used in the analysis. \code{gate} = "AND" or "OR", depending on whether a conjunctive o disjunctive nonparametric CDM is used. \code{pseudo.prob} = pseudo-posterior probability of each examinee mastering each attribute (experimental). \code{w.type} = weight type used for computing the pseudo-posterior probability (experimental): 1 = Power-of-2 weight; 2 = Exponential weight. \code{seed} = Numeric vector of length 1. NPS has a random component, so a seed is required for consistent results.
#' @param itemExposurecontrol  Scalar character. Item exposure control: \code{NULL} or progressive method (Barrada, Olea, Ponsoda, & Abad, 2008) with \code{"progressive"}. Default is \code{NULL}
#' @param b Scalar numeric. Acceleration parameter for the item exposure method. Only applies if \code{itemExposurecontrol = "progressive"}. In the progressive method the first item is selected at random and the last item (i.e., \code{MAXJ}) is selected purely based on \code{itemSelect}. The rest of the items are selected combining both a random and information components. The loss of importance of the random component will be linear with \code{b = 0}, inverse exponential with \code{b < 0}, or exponential with \code{b > 0}. Thus, \code{b} allows to optimize accuracy (\code{b < 0}) or item security (\code{b > 0}). Default is 2
#' @param maxr Scalar numeric. Value should be in the range 0-1. Maximum item exposure rate that is tolerated. Default is 1. Note that for \code{maxr < 1} parallel computing cannot be implemented
#' @param itemConstraint  Scalar character. Constraints that must be satisfied by the set of items applied: \code{NULL} or attribute constraint (Henson & Douglas, 2005) with \code{"attribute"}. If \code{"attribute"} is chosen, then each attribute must be measured at least a specific number of times indicated in the \code{constraint.args$ATTRIBUTEc} argument. Default is \code{NULL}
#' @param constraint.args A list of options when \code{itemConstraint != "NULL"}. At the moment it only includes the argument \code{ATTRIBUTEc} which must be a numeric vector of length \code{ncol(Q)} indicating the minimum number of items per attribute to be administered. Default is 3 
#' @param n.cores Scalar numeric. Number of cores to be used during parallelization. Default is 2
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
#' Robitzsch, A., Kiefer, T., George, A. C., & Uenlue, A. (2020). CDM: Cognitive Diagnosis Modeling. R package version 7.5-15. https://CRAN.R-project.org/package=CDM
#'
#' Yigit, H. D., Sorrel, M. A., de la Torre, J. (2018). Computerized adaptive testing for cognitively based multiple-choice data. \emph{Applied Psychological Measurement, 43}, 388-401.
#'
#' @import foreach
#' @import doSNOW
#' @import stats
#' @import utils
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
#'res.FIXJ.sum.real$recovery$plotPCV
#'res.FIXJ.sum.real$recovery$plotPCA
#'res.FIXJ.sum.real$item.exposure$plot
#'# VARJ summary
#'res.VARJ.sum.real <- cdcat.summary(cdcat.obj = res.VARJ, alpha = att)
#'res.VARJ.sum.real$recovery
#'res.VARJ.sum.real$CATlength$stats
#'res.VARJ.sum.real$CATlength$plot
#'res.VARJ.sum.real$item.exposure$plot
#'# vs. maximum observable accuracy
#'att.J <- GDINA::personparm(fit, "MAP")[, -(K+1)] # GDINA package
#'# att.J <- t(sapply(strsplit(as.character(fit$pattern$map.est), ""), as.numeric)) # CDM package
#'class.J <- GDINA::ClassRate(att, att.J) # upper-limit for accuracy
#'res.FIXJ.sum.obse <- cdcat.summary(cdcat.obj = res.FIXJ, alpha = att.J)
#'res.FIXJ.sum.obse$recovery$plotPCV + ggplot2::geom_hline(yintercept = class.J$PCV[K],
#'                                                         color = "firebrick3")
#'res.FIXJ.sum.obse$recovery$plotPCA + ggplot2::geom_hline(yintercept = class.J$PCA,
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
#'# small-scale assessment             #
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
#'                      NPS.args = list(Q = Q, gate = "AND", pseudo.prob = TRUE, w.type = 1,
#'                      seed = 12345))
#'res.NPS.VARJ <- cdcat(dat = dat, itemSelect = "NPS", FIXED.LENGTH = FALSE,
#'                      MAXJ = 25, precision.cut = 0.90, n.cores = 2,
#'                      NPS.args = list(Q = Q, gate = "AND", pseudo.prob = TRUE, w.type = 1,
#'                      seed = 12345))
#'
#'#---------------Results--------------#
#'res.NPS.FIXJ$est[[1]] # estimates for the first examinee (fixed-length)
#'res.NPS.VARJ$est[[1]] # estimates for the first examinee (fixed-precision)
#'att.plot(res.NPS.FIXJ, i = 1) # plot for estimates for the first examinee (fixed-length)
#'att.plot(res.NPS.VARJ, i = 1) # plot for estimates for the first examinee (fixed-precision)
#'# FIXJ summary
#'res.NPS.FIXJ.sum.real <- cdcat.summary(cdcat.obj = res.NPS.FIXJ, alpha = att) # vs. real accuracy
#'res.NPS.FIXJ.sum.real$recovery$plotPCV
#'res.NPS.FIXJ.sum.real$recovery$plotPCA
#'res.NPS.FIXJ.sum.real$item.exposure$plot
#'# VARJ summary
#'res.NPS.VARJ.sum.real <- cdcat.summary(cdcat.obj = res.NPS.VARJ, alpha = att)
#'res.NPS.VARJ.sum.real$recovery
#'res.NPS.VARJ.sum.real$CATlength$stats
#'res.NPS.VARJ.sum.real$CATlength$plot
#'res.NPS.VARJ.sum.real$item.exposure$plot
#'# vs. maximum observable accuracy
#'fit <- NPCD::AlphaNP(Y = dat, Q = Q, gate = "AND")
#'att.J <- fit$alpha.est
#'class.J <- GDINA::ClassRate(att, att.J) # upper-limit for accuracy
#'res.NPS.FIXJ.sum.obse <- cdcat.summary(cdcat.obj = res.NPS.FIXJ, alpha = att.J)
#'res.NPS.FIXJ.sum.obse$recovery$plotPCV + ggplot2::geom_hline(yintercept = class.J$PCV[K],
#'                                                             color = "firebrick3")
#'res.NPS.FIXJ.sum.obse$recovery$plotPCA + ggplot2::geom_hline(yintercept = class.J$PCA,
#'                                                             color = "firebrick3")
#'}
#' @export
#'
cdcat <- 
  function (fit = NULL, dat = NULL, itemSelect = "GDI", MAXJ = 20, FIXED.LENGTH = TRUE, 
            startRule = "random", startK = FALSE,
            att.prior = NULL, initial.distr = NULL, 
            precision.cut = 0.8, NPS.args = list(Q = NULL, gate = NULL, 
                                                 pseudo.prob = T, w.type = 1, seed = NULL),  
            itemExposurecontrol = NULL, b = 2, maxr = 1,
            itemConstraint = NULL, constraint.args = list(ATTRIBUTEc = NULL),
            n.cores = 2, print.progress = TRUE)
  {

    # Step 1: get data #############################################################################################
    Q <- NPS.args$Qg
    initial.distr.arg <- initial.distr
    if (!is.null(fit)) {
      if (!is.null(fit$extra$call)) {
        model <- ifelse(length(unique(fit$model)) == 1, unique(fit$model), 
                        "Combination")
        if (is.null(dat)) {
          X <- fit$options$dat
        }
        else {
          X <- dat
        }
        if (is.null(Q)) {
          Q <- fit$options$Q
        }
        LC.prob <- fit$LC.prob
      }
      else {
        model <- ifelse(length(unique(fit$rule)) == 1, unique(fit$rule), 
                        "Combination")
        if (is.null(dat)) {
          X <- fit$data
        }
        else {
          X <- dat
        }
        if (is.null(Q)) {
          Q <- fit$q.matrix
        }
        LC.prob <- matrix(NA, nrow = nrow(Q), ncol = 2^ncol(Q), 
                          dimnames = list(paste("Item", 1:nrow(Q)), 
                                          apply(GDINA::attributepattern(ncol(Q)), 1, 
                                                paste, collapse = "")))
        for (j in 1:nrow(Q)) {
          kj <- which(Q[j, ] == 1)
          prob.j.l <- fit$probitem[fit$probitem$itemno == 
                                     j, ]$prob
          tmp <- factor(apply(GDINA::attributepattern(ncol(Q))[, 
                                                               kj, drop = F], 1, paste, collapse = ""))
          tmp <- factor(tmp, levels = apply(unique(GDINA::attributepattern(length(kj))), 
                                            1, paste, collapse = ""))
          for (l in 1:length(prob.j.l)) {
            LC.prob[j, which(as.numeric(tmp) == l)] <- prob.j.l[l]
          }
        }
      }
      LC.prob[LC.prob > 1] <- 1
      LC.prob[LC.prob < 0] <- 0
    } else {
      model <- paste0("NP_", NPS.args$gate)
      Q <- NPS.args$Q
      X <- dat
    }
    N <- nrow(X)
    q.matrix <- Q
    J <- nrow(Q)
    K <- ncol(Q)
    L <- 2^K
    pattern <- GDINA::attributepattern(K)
    Lclass <- apply(GDINA::attributepattern(K), MARGIN = 1, FUN = function(x) {
      paste(x, collapse = "")
    })
    if (is.null(att.prior)) {
      att.prior <- rep(1/L, L)
    }
    if (is.null(initial.distr.arg)) {
      initial.distr <- matrix(data = rep(rep(1/L, L), N), nrow = N)
    }
    if (!is.null(initial.distr.arg)) {
      initial.distr <- matrix(data = rep(initial.distr, N), 
                              nrow = N, byrow = TRUE)
    }
    
    # Step 2: check errors and warnings #############################################################################################
    if (!(class(fit) != "GDINA" | class(fit) != "gdina")) {
      stop("fit must be an object of class 'GDINA' or 'gdina'")
    }
    if (!(itemSelect %in% c("GDI", "JSD", "MPWKL", 
                            "PWKL", "NPS", "random"))) {
      stop("itemSelect must be one of the following: 'GDI', 'JSD', 'MPWKL', 'PWKL', 'NPS', 'random'")
    }
    if ((maxr > 1) | (maxr < 1)) {
      stop("maxr should be in the range 0-1")
    }
    if (is.null(itemConstraint) == FALSE) {
      if (is.null(constraint.args$ATTRIBUTEc)) {
        constraint.args$ATTRIBUTEc <- 3}
    }
    if (itemSelect == "NPS") {
      if (maxr != 1) {
        warning("maxr argument does not operate  with itemSelect == 'NPS'")
      }
      if (!is.null(itemExposurecontrol)) {
        warning("itemExposurecontrol is set to 'NULL' with itemSelect == 'NPS'")
      }
      if (is.null(dat)) {
        stop("dat required when itemSelect == 'NPS'")
      }
      if (is.null(Q)) {
        stop("Q required when itemSelect == 'NPS'")
      }
      if (is.null(NPS.args$pseudo.prob)) {
        NPS.args$pseudo.prob <- TRUE
      }
      if (is.null(NPS.args$w.type)) {
        NPS.args$w.type <- 1
      }
      if (is.null(NPS.args$seed)) {
        NPS.args$seed <- sample(1:1e+06, size = 1)
      }
      if (is.null(NPS.args$gate)) {
        stop("gate == 'AND' or gate == 'OR' required in NPS.args")
      }     else {
        if (!(NPS.args$gate %in% c("AND", "OR"))) {
          stop("gate == 'AND' or gate == 'OR' required in NPS.args")
        }
      }
      if (!(NPS.args$pseudo.prob %in% c(TRUE, FALSE))) {
        stop("pseudo.prob == TRUE or pseudo.prob == FALSE required in NPS.args")
      }
      if (!FIXED.LENGTH & !NPS.args$pseudo.prob) {
        warning("FIXED.LENGHT = FALSE is not available with itemSelect == 'NPS' if NPS.args$pseudo.prob == FALSE. NPS.args$pseudo.prob = TRUE applied instead")
      }
    }   else {
      if (is.null(fit)) {
        stop("fit required when itemSelect != 'NPS'")
      }
    }
    if (n.cores > parallel::detectCores()) {
      stop("n.cores cannot be higher than the number of cores in the CPU")
    }
    if (length(att.prior) != 2^K) {
      stop("att.prior should be of length 2^K")
    }
    if (sum(att.prior) != 1) {
      stop("att.prior should add up to 1")
    }
    if(!is.null(itemExposurecontrol)) {
      if (itemExposurecontrol == "restricprogressive") {
        if (maxr != 1) {
          stop("Restrictive progressive item exposure control cannot be implemented with maxr != 1. Set maxr = 1")
        }
      } 
    }
    
    # Step 3: conduct CD-CAT #############################################################################################
    out <- list()
    item.usage <- NULL
    
    if (itemSelect != "NPS") {
      
      if(maxr == 1) {
        
        cl <- parallel::makeCluster(n.cores, type = "SOCK")
        doSNOW::registerDoSNOW(cl)
        
        if (print.progress == TRUE) {
          pb <- txtProgressBar(max = N, style = 3)
          progress <- function(n) setTxtProgressBar(pb, n)
          opts <- list(progress = progress)
        }
        
        if (print.progress == FALSE) {
          opts <- NULL
        }
        
        out <- foreach(i = 1:N, .options.snow = opts,
                       .export = c("GDI.M", "proggresive.f", "H", "JSD.DICO.M", "PWKL.M", "MPWKL.M"),
                       .inorder = T) %dopar% {
                         try({
                           
                           # initialize
                           est.cat <- matrix(NA, nrow = 1, ncol = 10 + K, 
                                             dimnames = list(1, c("j", "qj", "xj",
                                                                  "ML", "nmodesML", "Lik", 
                                                                  "MAP", "nmodesMAP", "Post", 
                                                                  "EAP", paste("K", 1:K, sep = ""))))[-1, 
                                                                  ]
                           item.log <- NULL
                           
                           # get data
                           item.usage <- lapply(out, function(x) x[[2]])
                           mlogPost_GDI <- initial.distr[i, ]
                           
                           GDI <- switch(itemSelect, 
                                         GDI = GDI.M(LC.prob = LC.prob, mlogPost_GDI), 
                                         JSD = JSD.DICO.M(LC.prob = LC.prob, mlogPost_GDI), 
                                         MPWKL = MPWKL.M(LC.prob = LC.prob, mlogPost_GDI), 
                                         PWKL = PWKL.M(LC.prob = LC.prob, mlogPost_GDI, point.est = sample(which(mlogPost_GDI == max(mlogPost_GDI)), 1)), 
                                         random = runif(J, 0, 1)) 
                           
                           if(!is.null(itemExposurecontrol)) {
                             if (itemExposurecontrol == "progressive") {
                               GDI <- proggresive.f(jjselect = 1, MAXJ = MAXJ, b = b, GDI = GDI)
                             } 
                           }
                           names(GDI) <- 1:length(GDI)
                           
                           # CD-CAT
                           jjselect <- 1
                           jjcatGDI <- IP <- NULL
                           
                           for (jjselect in 1:MAXJ) {
                             
                             if(startK){
                               if(jjselect == 1){
                                 startK.pool <- which(rowSums(Q) == 1)
                                 if(startRule == "random"){
                                   jjcatGDI[jjselect] <- as.numeric(sample(startK.pool, 1))
                                 }
                                 if(startRule == "max"){
                                   jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI[startK.pool])))
                                 }
                               } else if(jjselect <= K){
                                 past.k <- which(colSums(Q[jjcatGDI,,drop = FALSE]) > 0)
                                 startK.pool <- which(rowSums(Q[, past.k, drop = FALSE]) == 0 & rowSums(Q[, -past.k, drop = FALSE]) == 1)
                                 jjcatGDI[jjselect] <- as.numeric(sample(startK.pool, 1)) 
                               } else if(jjselect > K){
                                 jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI[eligible])))
                               }
                             } else {
                               if(jjselect == 1){
                                 if(startRule == "random"){
                                   jjcatGDI[jjselect] <- as.numeric(sample(1:J, 1))
                                 }
                                 if(startRule == "max"){
                                   jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI)))
                                 }
                               } else {
                                 jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI[eligible])))
                               }
                             }
                             
                             x.jj <- X[i, jjcatGDI[jjselect]]
                             IP <- cbind(IP, x.jj)
                             p.xi <- matrix(LC.prob[jjcatGDI, ], ncol = L)
                             p.xi.aj1_GDI <- matrix(0, nrow = L, ncol = length(IP))
                             
                             for (l in 1L:L) {
                               for (jj in 1L:length(IP)) {
                                 p.xi.aj1_GDI[l, jj] <- (p.xi[jj, l]^IP[jj]) * ((1 - p.xi[jj, l])^(1 - IP[jj]))
                               }
                             }
                             
                             p.xi.aj_GDI <- (t(as.matrix(apply(p.xi.aj1_GDI, 1, prod))))
                             mPost_GDI <- p.xi.aj_GDI * att.prior
                             msumPost_GDI <- sum(mPost_GDI)
                             mlogPost_GDI <- mPost_GDI/msumPost_GDI
                             p.aj_GDI <- mlogPost_GDI %*% pattern
                             
                             # check constraints
                             if(!is.null(itemConstraint)) {
                               # Attribute constraint. Note: eligible items are those who have not being applied and measure an attribute that has not reach a limit
                               Q.tmp <- as.matrix(q.matrix[jjcatGDI, ])
                               if(nrow(Q.tmp) == 1) {Q.tmp <- matrix(Q.tmp, nrow = 1)}
                               Q.tmp.AC <- colSums(Q.tmp) < constraint.args$ATTRIBUTEc # weight for item selection
                               if (sum(Q.tmp.AC) == 1) {
                                 cond2 <- (Q[, which(Q.tmp.AC == TRUE)] > 0) 
                               } else {
                                 cond2 <- (rowSums(Q[, which(Q.tmp.AC == TRUE)]) > 0)
                               }
                               if (sum((!(1:nrow(Q)%in%jjcatGDI)) & cond2) == 0) {
                                 eligible <- which((!(1:nrow(Q)%in%jjcatGDI)))
                               } else {
                                 eligible <- (1:nrow(Q))[(!(1:nrow(Q)%in%jjcatGDI)) & cond2] 
                               }
                             } else {
                               eligible <- which((!(1:nrow(Q)%in%jjcatGDI)))
                             }
                             
                             # update itemSelect
                             GDI <- switch(itemSelect, 
                                           GDI = GDI.M(LC.prob = LC.prob, mlogPost_GDI), 
                                           JSD = JSD.DICO.M(LC.prob = LC.prob, mlogPost_GDI), 
                                           MPWKL = MPWKL.M(LC.prob = LC.prob, mlogPost_GDI), 
                                           PWKL = PWKL.M(LC.prob = LC.prob, mlogPost_GDI, point.est = sample(which(mlogPost_GDI == max(mlogPost_GDI)), 1)), 
                                           random = runif(J, 0, 1)) 
                             if(!is.null(itemExposurecontrol)) {
                               if (itemExposurecontrol == "progressive") {
                                 GDI <- proggresive.f(jjselect = jjselect + 1, MAXJ = MAXJ, b = b, GDI = GDI)
                               } 
                             }
                             names(GDI) <- 1:length(GDI)
                             
                             # output
                             est.cat <- rbind(est.cat, c(jjcatGDI[jjselect], 
                                                         paste(q.matrix[jjcatGDI[jjselect], ], collapse = ""), 
                                                         IP[jjselect],
                                                         Lclass[which.max(p.xi.aj_GDI)], length(which(p.xi.aj_GDI == max(p.xi.aj_GDI))), 
                                                         round(max(p.xi.aj_GDI), 5), 
                                                         Lclass[which.max(mlogPost_GDI)], length(which(mlogPost_GDI == max(mlogPost_GDI))), 
                                                         round(max(mlogPost_GDI), 5), 
                                                         paste(as.numeric(p.aj_GDI > 0.5), collapse = ""), 
                                                         round(p.aj_GDI, 5)))
                             
                             # check stopping rule
                             jjselect <- jjselect + 1
                             if ((jjselect) > MAXJ) {
                               break
                             }
                             if (!FIXED.LENGTH) {
                               if (max(mlogPost_GDI) >= precision.cut) {
                                 break
                               }
                             }
                           }
                           
                           item.log <- c(item.log, jjcatGDI)
                           est.cat <- as.data.frame(est.cat)
                           out[[i]] <- list(est.cat = as.data.frame(est.cat), 
                                            item.usage = item.log)
                         })
                       }
        
      } else { # end maxr == 1
        
        for (i in 1:N){
          
          if(isTRUE(print.progress)) {  cat('\r Processing examinee', i, 'of', N)}
          
          try({
            # initialize
            est.cat <- matrix(NA, nrow = 1, ncol = 10 + K, 
                              dimnames = list(1, c("j", "qj", "xj",
                                                   "ML", "nmodesML", "Lik", 
                                                   "MAP", "nmodesMAP", "Post", 
                                                   "EAP", paste("K", 1:K, sep = ""))))[-1, 
                                                   ]
            item.log <- NULL
            
            # get data
            item.usage <- lapply(out, function(x) x[[2]])
            mlogPost_GDI <- initial.distr[i, ]
            
            GDI <- switch(itemSelect, 
                          GDI = GDI.M(LC.prob = LC.prob, mlogPost_GDI), 
                          JSD = JSD.DICO.M(LC.prob = LC.prob, mlogPost_GDI), 
                          MPWKL = MPWKL.M(LC.prob = LC.prob, mlogPost_GDI), 
                          PWKL = PWKL.M(LC.prob = LC.prob, mlogPost_GDI, point.est = sample(which(mlogPost_GDI == max(mlogPost_GDI)), 1)), 
                          random = runif(J, 0, 1)) 
            if(!is.null(itemExposurecontrol)) {
              if (itemExposurecontrol == "progressive") {
                GDI <- proggresive.f(jjselect = 1, MAXJ = MAXJ, b = b, GDI = GDI)
                }
              }
            
            # CD-CAT
            jjselect <- 1
            jjcatGDI <- IP <- NULL
            
            usage <-c(unlist(item.usage),jjcatGDI)
            r_i <- as.vector(table(factor(usage, levels = 1: nrow(Q))) / N)
            tmp <- rep(1, times = length(GDI))
            tmp[which((1 - r_i/maxr) < 0)] <- 0
            GDI <- tmp * GDI  # note: items with r_i > maxr cannot be applied anymore
            names(GDI) <- 1:length(GDI)
            
            for (jjselect in 1:MAXJ) {
              
              if(startK){
                if(jjselect == 1){
                  startK.pool <- which(rowSums(Q) == 1)
                  if(startRule == "random"){
                    jjcatGDI[jjselect] <- as.numeric(sample(startK.pool, 1))
                  }
                  if(startRule == "max"){
                    jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI[startK.pool])))
                  }
                } else if(jjselect <= K){
                  past.k <- which(colSums(Q[jjcatGDI,,drop = FALSE]) > 0)
                  startK.pool <- which(rowSums(Q[, past.k, drop = FALSE]) == 0 & rowSums(Q[, -past.k, drop = FALSE]) == 1)
                  jjcatGDI[jjselect] <- as.numeric(sample(startK.pool, 1)) 
                } else if(jjselect > K){
                  jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI[eligible])))
                }
              } else {
                if(jjselect == 1){
                  if(startRule == "random"){
                    jjcatGDI[jjselect] <- as.numeric(sample(1:J, 1))
                  }
                  if(startRule == "max"){
                    jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI)))
                  }
                } else {
                  jjcatGDI[jjselect] <- as.numeric(names(which.max(GDI[eligible])))
                }
              }
              
              x.jj <- X[i, jjcatGDI[jjselect]]
              IP <- cbind(IP, x.jj)
              p.xi <- matrix(LC.prob[jjcatGDI, ], ncol = L)
              p.xi.aj1_GDI <- matrix(0, nrow = L, ncol = length(IP))
              for (l in 1L:L) {
                for (jj in 1L:length(IP)) {
                  p.xi.aj1_GDI[l, jj] <- (p.xi[jj, l]^IP[jj]) * ((1 - p.xi[jj, l])^(1 - IP[jj]))
                }
              }
              p.xi.aj_GDI <- (t(as.matrix(apply(p.xi.aj1_GDI, 1, prod))))
              mPost_GDI <- p.xi.aj_GDI * att.prior
              msumPost_GDI <- sum(mPost_GDI)
              mlogPost_GDI <- mPost_GDI/msumPost_GDI
              p.aj_GDI <- mlogPost_GDI %*% pattern
              
              # check constraints
              if(is.null(itemConstraint) == FALSE) {
                # Attribute constraint. Note: eligible items are those who have not being applied and measure an attribute that has not reach a limit
                Q.tmp <- as.matrix(q.matrix[jjcatGDI, ])
                if(nrow(Q.tmp) == 1) {Q.tmp <- matrix(Q.tmp, nrow = 1)}
                Q.tmp.AC <- colSums(Q.tmp) < constraint.args$ATTRIBUTEc # weight for item selection
                if (sum(Q.tmp.AC) == 1) {
                  cond2 <- (Q[, which(Q.tmp.AC == TRUE)] > 0) 
                } else {
                  cond2 <- (rowSums(Q[, which(Q.tmp.AC == TRUE)]) > 0)
                }
                if (sum((!(1:nrow(Q)%in%jjcatGDI)) & cond2) == 0) {
                  eligible <- which((!(1:nrow(Q)%in%jjcatGDI)))
                } else {
                  eligible <- (1:nrow(Q))[(!(1:nrow(Q)%in%jjcatGDI)) & cond2] 
                }
              } else {
                eligible <- which((!(1:nrow(Q)%in%jjcatGDI)))
              }
              
              # update itemSelect
              GDI <- switch(itemSelect, 
                            GDI = GDI.M(LC.prob = LC.prob, mlogPost_GDI), 
                            JSD = JSD.DICO.M(LC.prob = LC.prob, mlogPost_GDI), 
                            MPWKL = MPWKL.M(LC.prob = LC.prob, mlogPost_GDI), 
                            PWKL = PWKL.M(LC.prob = LC.prob, mlogPost_GDI, point.est = sample(which(mlogPost_GDI == max(mlogPost_GDI)), 1)), 
                            random = runif(J, 0, 1)) 
              if(!is.null(itemExposurecontrol)) {
                if (itemExposurecontrol == "progressive") {
                  GDI <- proggresive.f(jjselect = jjselect + 1, MAXJ = MAXJ, b = b, GDI = GDI)
                } 
              }
              usage <-c(unlist(item.usage),jjcatGDI)
              r_i <- as.vector(table(factor(usage, levels = 1: nrow(Q))) / N)
              tmp <- rep(1, times = length(GDI))
              tmp[which((1 - r_i/maxr) < 0)] <- 0
              GDI <- tmp * GDI  
              names(GDI) <- 1:length(GDI)
              
              # output
              est.cat <- rbind(est.cat, c(jjcatGDI[jjselect], 
                                          paste(q.matrix[jjcatGDI[jjselect], ], collapse = ""), 
                                          IP[jjselect],
                                          Lclass[which.max(p.xi.aj_GDI)], length(which(p.xi.aj_GDI == max(p.xi.aj_GDI))), 
                                          round(max(p.xi.aj_GDI), 5), 
                                          Lclass[which.max(mlogPost_GDI)], length(which(mlogPost_GDI == max(mlogPost_GDI))), 
                                          round(max(mlogPost_GDI),                                                                                                                                                                                                   5), paste(as.numeric(p.aj_GDI > 0.5), collapse = ""), 
                                          round(p.aj_GDI, 5)))
              
              # check stopping rule
              jjselect <- jjselect + 1
              if ((jjselect) > MAXJ) {
                break
              }
              if (!FIXED.LENGTH) {
                if (max(mlogPost_GDI) >= precision.cut) {
                  break
                }
              }
            }
            item.log <- c(item.log, jjcatGDI)
            est.cat <- as.data.frame(est.cat)
            out[[i]] <- list(est.cat = as.data.frame(est.cat), 
                             item.usage = item.log)
          })
        }
        
      } # end maxr != 1
      
    } else { # end itemSelect != "NPS"
      
        cl <- parallel::makeCluster(n.cores, type = "SOCK")
        doSNOW::registerDoSNOW(cl)
        
        if (print.progress == TRUE) {
          pb <- txtProgressBar(max = N, style = 3)
          progress <- function(n) setTxtProgressBar(pb, n)
          opts <- list(progress = progress)
        }
        
        if (print.progress == FALSE) {
          opts <- NULL
        }

      # initialize
      att.pat <- matrix(NA, nrow = N, ncol = K)
      exp.items <- list()
      resp.pat <- list()
      
      # get data
      if (NPS.args$gate == "AND") {
        gate <- rep("AND", J)
      }
      if (NPS.args$gate == "OR") {
        gate <- rep("OR", J)
      }
      if (!FIXED.LENGTH & !NPS.args$pseudo.prob) {
        NPS.args$pseudo.prob <- TRUE
      }
      
      # CD-CAT
      out <- foreach(i = 1:N, .options.snow = opts, .packages = "NPCD", 
                     .export = c("pseudoP", "NPC.eta"), .inorder = T) %dopar% {
                       try({
                         if (!is.null(NPS.args$seed)) {
                           set.seed(NPS.args$seed + i)
                         }
                         if (NPS.args$pseudo.prob) {
                           out.i <- matrix(NA, ncol = 8 + K)[-1, ]
                           colnames(out.i) <- c("j", "qj", 
                                                "xj", "alpha", "HD.alpha", 
                                                "alpha2", "HD.alpha2", "HD.diff", 
                                                paste0("pP.K", 1:K))
                         }
                         else {
                           out.i <- matrix(NA, ncol = 8)[-1, ]
                           colnames(out.i) <- c("j", "qj", 
                                                "xj", "alpha", "HD.alpha", 
                                                "alpha2", "HD.alpha2", "HD.diff")
                         }
                         x <- X[i, ]
                         y <- c()
                         Qy <- c()
                         gatey <- c()
                         R0 <- as.data.frame(Q)
                         gate0 <- gate
                         e1 <- c(1, rep(0, K - 1))
                         hit <- which(sapply(apply(R0, 1, function(x) which(all(as.double(x) == 
                                                                                  e1))), sum) == 1)
                         if (length(hit) > 1) {
                           hit <- sample(hit, 1)
                         }
                         R1 <- R0[-hit, ]
                         gate1 <- gate0[-hit]
                         y <- c(y, x[hit])
                         x <- x[-hit]
                         Qy <- rbind(Qy, R0[hit, ])
                         gatey <- c(gatey, gate0[hit])
                         if (NPS.args$pseudo.prob) {
                           out.i <- rbind(out.i, c(as.double(rownames(Qy)), 
                                                   paste(Qy[1, ], collapse = ""), y, 
                                                   rep(NA, 5 + K)))
                         } else {
                           out.i <- rbind(out.i, c(as.double(rownames(Qy)), 
                                                   paste(Qy[1, ], collapse = ""), y, 
                                                   rep(NA, 5)))
                         }
                         
                         for (k in 2:K) {
                           hit <- NULL
                           assign(paste0("e", k), c(rep(0, k - 
                                                          1), 1, rep(0, K - k)))
                           while (length(hit) == 0) {
                             qj.prev <- matrix(get(paste0("e", 
                                                          k)), ncol = K, byrow = T)
                             for (w in 1:(k - 1)) {
                               qj.prev <- rbind(qj.prev, sample(c(0, 
                                                                  1), 1) * get(paste0("e", w)) * 
                                                  y[w])
                             }
                             qj <- colSums(qj.prev)
                             qj[qj > 0] <- 1
                             hit <- as.vector(which(sapply(apply(get(paste0("R", 
                                                                            k - 1)), 1, function(x) which(all(as.double(x) == 
                                                                                                                qj))), sum) == 1))
                             if (length(hit) > 1) {
                               hit <- sample(hit, 1)
                             }
                           }
                           assign(paste0("R", k), get(paste0("R", k - 1))[-hit, ])
                           assign(paste0("gate", k), get(paste0("gate", k - 1))[-hit])
                           y <- c(y, x[hit])
                           x <- x[-hit]
                           Qy <- rbind(Qy, get(paste0("R", k - 1))[hit, ])
                           gatey <- c(gatey, get(paste0("gate", k - 1))[hit])
                           if (NPS.args$pseudo.prob) {
                             out.i <- rbind(out.i, c(as.double(rownames(Qy))[k], 
                                                     paste(Qy[k, ], collapse = ""), 
                                                     y[k], rep(NA, 5 + K)))
                           } else {
                             out.i <- rbind(out.i, c(as.double(rownames(Qy))[k], 
                                                     paste(Qy[k, ], collapse = ""), 
                                                     y[k], rep(NA, 5)))
                           }
                         }
                         alphanp <- NPCD::AlphaNP(t(as.matrix(y)), Qy, gatey[1])
                         pat.dist <- cbind(alphanp$pattern, alphanp$loss.matrix)
                         if (K > 1) {
                           random.order <- sample(1:K, K, F)
                         }
                         for (k in random.order) {
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
                         if (NPS.args$pseudo.prob) {
                           pP.k <- as.vector(pseudoP(alphanp, inv.type = "min", w.type = NPS.args$w.type))
                           out.i[K, 9:(8 + K)] <- round(pP.k, 5)
                         }
                         t <- K + 1
                         while (t <= MAXJ) {
                           hit <- NULL
                           while (is.null(hit)) {
                             items <- sample(1:nrow(get(paste0("R", t - 1))))
                             for (j in items) {
                               qj <- as.double(get(paste0("R", t - 1))[j, ])
                               eta.a <- NPC.eta(qj, get(paste0("a.", t - 1)), 
                                                get(paste0("gate", t - 1))[j])
                               eta.a2 <- NPC.eta(qj, get(paste0("a2.", t - 1)), 
                                                 get(paste0("gate", t - 1))[j])
                               if (eta.a != eta.a2) {
                                 hit <- as.vector(which(sapply(apply(get(paste0("R", t - 1)), 1, 
                                                                     function(x) which(all(as.double(x) == qj))), sum) == 1))
                                 if (length(hit) > 1) {
                                   hit <- sample(hit, 1)
                                 }
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
                             if (is.null(hit)) {
                               m <- m + 1
                               alphanp <- NPCD::AlphaNP(t(as.matrix(y)), Qy, gatey[1])
                               pat.dist <- cbind(alphanp$pattern, alphanp$loss.matrix)
                               if (K > 1) {
                                 random.order <- sample(1:K, K, F)
                               }
                               for (k in random.order) {
                                 pat.dist <- pat.dist[order(pat.dist[, k]), ]
                               }
                               pat.dist <- pat.dist[order(pat.dist[, K + 1]), ]
                               assign(paste0("a2.", t - 1), pat.dist[order(pat.dist[, K + 1]), ][m, 1:K])
                             }
                           }
                           alphanp <- NPCD::AlphaNP(t(as.matrix(y)), Qy, gatey[1])
                           pat.dist <- cbind(alphanp$pattern, alphanp$loss.matrix)
                           if (K > 1) {
                             random.order <- sample(1:K, K, F)
                           }
                           for (k in random.order) {
                             pat.dist <- pat.dist[order(pat.dist[, k]), ]
                           }
                           pat.dist <- pat.dist[order(pat.dist[, K + 1]), ]
                           assign(paste0("a.", t), pat.dist[order(pat.dist[, K + 1]), ][1, 1:K])
                           m <- 2
                           assign(paste0("a2.", t), pat.dist[order(pat.dist[, K + 1]), ][m, 1:K])
                           if (!NPS.args$pseudo.prob) {
                             out.i <- rbind(out.i, c(as.double(rownames(Qy))[t], 
                                                     paste(Qy[t, ], collapse = ""), 
                                                     y[t], paste(get(paste0("a.", t)), collapse = ""), 
                                                     pat.dist[order(pat.dist[, K + 1]), ][1, K + 1], 
                                                     paste(get(paste0("a2.", t)), collapse = ""), 
                                                     pat.dist[order(pat.dist[, K + 1]), ][2, K + 1], 
                                                     pat.dist[order(pat.dist[, K + 1]), ][2, K + 1] - pat.dist[order(pat.dist[, K + 1]), ][1, K + 1]))
                           } else {
                             if (NPS.args$pseudo.prob) {
                               pP.k <- as.vector(pseudoP(alphanp, inv.type = "min", w.type = NPS.args$w.type))
                               out.i <- rbind(out.i, c(as.double(rownames(Qy))[t], 
                                                       paste(Qy[t, ], collapse = ""), 
                                                       y[t], paste(get(paste0("a.", t)), collapse = ""), 
                                                       pat.dist[order(pat.dist[, K + 1]), ][1, K + 1], 
                                                       paste(get(paste0("a2.", t)), collapse = ""), 
                                                       pat.dist[order(pat.dist[, K + 1]), ][2, K + 1], 
                                                       pat.dist[order(pat.dist[, K + 1]), ][2, K + 1] - pat.dist[order(pat.dist[, K + 1]), ][1, K + 1], 
                                                       round(pP.k, 5)))
                             }
                           }
                           if (NPS.args$pseudo.prob) {
                             if (!FIXED.LENGTH) {
                               cut.pP.k <- pP.k
                               cut.pP.k[cut.pP.k < 0.5] <- 1 - cut.pP.k[cut.pP.k < 0.5]
                               if (all(cut.pP.k > precision.cut)) {
                                 break
                               }
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
                         if (NPS.args$pseudo.prob) {
                           for (k in 1:K) {
                             out.i[, 8 + k] <- as.double(as.character(out.i[, 8 + k]))
                           }
                         }
                         out[[i]] <- list(est.cat = out.i, item.usage = exp.items[[i]])
                       })
                     }
    } # end itemSelect == "NPS"
    
    if ((print.progress == TRUE) & (maxr == 1)) {
      close(pb)
      parallel::stopCluster(cl)
    }
    
    res <- list()
    res$est <- out
    res$specifications <- list(fit = fit, dat = dat, Q = Q, model = model, itemSelect = itemSelect, MAXJ = MAXJ, FIXED.LENGTH = FIXED.LENGTH, 
                               startK = startK, startRule = startRule, att.prior = att.prior, initial.distr = initial.distr.arg, 
                               precision.cut = precision.cut, NPS.args = NPS.args, 
                               itemExposurecontrol = itemExposurecontrol, b = b, maxr = maxr,
                               itemConstraint = itemConstraint, constraint.args = constraint.args,
                               n.cores = n.cores, print.progress = TRUE)
    class(res) <- "cdcat"
    return(res)
  }