#' Data generation
#'
#' This function can be used to generate datasets based on object of class \code{gen.itembank}.
#' The user can manipulate the examinees' attribute distribution or provide a matrix of attribute profiles.
#' Data are simulated using the \code{GDINA::simGDINA} function (Ma & de la Torre, 2020).
#'
#' @param N Numeric vector of length 1. Sample size for the datasets
#' @param R Numeric vector of length 1. Number of datasets replications. Default is 1
#' @param item.bank An object of class \code{gen.itembank}
#' @param att.profiles Numeric matrix indicating the true attribute profile for each examinee (\emph{N} examinees x \emph{K} attributes). If \code{NULL} (by default), \code{att.dist} must be specified
#' @param att.dist Numeric vector of length 2^\emph{K}, where \emph{K} is the number of attributes. Distribution for attribute simulation. It can be \code{"uniform"} (by default), \code{"higher.order"}, \code{"mvnorm"}, or \code{"categorical"}. See \code{simGDINA} function of package \code{GDINA} for more information. Only used when \code{att.profiles = NULL}
#' @param mvnorm.parm A list of arguments for multivariate normal attribute distribution (\code{att.dist = "mvnorm"}). See \code{simGDINA} function of package \code{GDINA} for more information
#' @param higher.order.parm A list of arguments for higher-order attribute distribution (\code{att.dist = "mvnorm"}). See \code{simGDINA} function of package \code{GDINA} for more information
#' @param categorical.parm A list of arguments for categorical attribute distribution (\code{att.dist = "mvnorm"}). See \code{simGDINA} function of package \code{GDINA} for more information
#' @param seed Numeric vector of length 1. A scalar to use with \code{set.seed}
#'
#' @return \code{gen.itembank} returns an object of class \code{gen.itembank}.
#' \describe{
#' \item{simdat}{An array containing the simulated responses (dimensions N examinees x J items x R replicates). If \code{R = 1}, a matrix is provided}
#' \item{simalpha}{An array containing the simulated attribute profiles (dimensions N examinees x K attributes x R replicates). If \code{R = 1}, a matrix is provided}
#' \item{specifications}{A list that contains all the specifications}
#' }
#'
#' @references
#'
#' Ma, W. & de la Torre, J. (2020). GDINA: The generalized DINA model framework. R package version 2.7.9. Retrived from https://CRAN.R-project.org/package=GDINA
#'
#' @examples
#'####################################
#'# Example 1.                       #
#'# Generate dataset (GDINA item     #
#'# parameters and uniform attribute #
#'# distribution)                    #
#'####################################
#'
#' Q <- sim180GDINA$simQ
#' bank <- gen.itembank(Q = Q, mean.IQ = .70, range.IQ = .20, model = "GDINA")
#'
#' simdata <- gen.data(N = 1000, item.bank = bank)
#'
#'####################################
#'# Example 2.                       #
#'# Generate multiple datasets (DINA #
#'# model and multivariate normal    #
#'# attribute distribution)          #
#'####################################
#'
#'Q <- sim180GDINA$simQ
#'K <- ncol(Q)
#'bank <- gen.itembank(Q = Q, mean.IQ = .70, range.IQ = .20, model = "DINA")
#'
#'cutoffs <- qnorm(c(1:K)/(K+1))
#'m <- rep(0,K)
#'vcov <- matrix(0.5,K,K)
#'diag(vcov) <- 1
#'simdata <- gen.data(N = 1000, R = 20, item.bank = bank, att.dist = "mvnorm",
#'                    mvnorm.parm = list(mean = m, sigma = vcov, cutoffs = cutoffs))
#'
#'####################################
#'# Example 3.                       #
#'# Generate dataset (multiple       #
#'# models and higher-order          #
#'# attribute distribution)          #
#'####################################
#'
#'Q <- sim180GDINA$simQ
#'K <- ncol(Q)
#'model <- sample(c("DINA", "DINO", "ACDM"), size = nrow(Q), replace = TRUE)
#'bank <- gen.itembank(Q = Q, mean.IQ = .70, range.IQ = .20, model = model)
#'
#'N <- 1000
#'theta <- rnorm(N)
#'lambda <- data.frame(a = runif(K, 0.7, 1.3), b = seq( -2, 2, length.out = K))
#'simdata <- gen.data(N = N, item.bank = bank, att.dist = "higher.order",
#'                    higher.order.parm = list(theta = theta,lambda = lambda))
#'
#'####################################
#'# Example 4.                       #
#'# Generate dataset (GDINA model    #
#'# and given attribute profiles)    #
#'####################################
#'
#'Q <- sim180GDINA$simQ
#'K <- ncol(Q)
#'bank <- gen.itembank(Q = Q, mean.IQ = .70, range.IQ = .20, model = "GDINA")
#'
#'att.profiles <- matrix(data = c(1,0,0,0,0,
#'                                1,1,0,0,0,
#'                                1,1,1,0,0,
#'                                1,1,1,1,1), ncol = K, byrow = TRUE)
#'simdata <- gen.data(item.bank = bank, att.profiles = att.profiles)
#'
#' @export
#'
gen.data <- function(N = NULL, R = 1, item.bank = NULL, att.profiles = NULL, att.dist = "uniform",
                     mvnorm.parm = list(mean = NULL, sigma = NULL, cutoffs = NULL),
                     higher.order.parm = list(theta = NULL, lambda = NULL),
                     categorical.parm = list(att.prior = NULL),
                     seed = NULL)
  {

  #----------------------------
  # Warning end error messages
  #----------------------------

  if(!is.null(att.profiles)){N <- nrow(att.profiles)}
  if(is.null(N)){stop("N argument missing, with no default.")}
  if(R < 1){stop("R argument must be a positive scalar.")}
  if(is.null(item.bank)){stop("item.bank argument missing, with no default.")}
  if(!(att.dist %in% c("uniform", "mvnorm", "higher.order", "categorical"))){stop("att.dist must be equal to 'uniform', 'mvnorm', 'higher.order', or 'categorical'.")}
  if(att.dist == "mvnorm"){
    if(is.null(mvnorm.parm$mean)){stop("mvnorm.parm$mean argument missing, with no default.")}
    if(is.null(mvnorm.parm$sigma)){stop("mvnorm.parm$sigma argument missing, with no default.")}
    if(is.null(mvnorm.parm$cutoffs)){stop("mvnorm.parm$cutoffs argument missing, with no default.")}
  } else if(att.dist == "higher.order"){
    if(is.null(higher.order.parm$theta)){stop("higher.order.parm$theta argument missing, with no default.")}
    if(is.null(higher.order.parm$lambda)){stop("higher.order.parm$lambda argument missing, with no default.")}
  } else if(att.dist == "categorical"){
    if(is.null(categorical.parm$att.prior)){stop("categorical.parm$att.prior argument missing, with no default.")}
  }

  #----------------
  # Gather objects
  #----------------

  if(!is.null(seed)){set.seed(seed)}
  if(is.null(seed)){seed <- sample(1:1000000, size = 1)}
  if(is.null(item.bank$simQ)){
    Q <- item.bank$specifications$Q
  } else {
    Q <- item.bank$simQ
  }
  J <- nrow(Q)
  K <- ncol(Q)
  catprob.parm <- item.bank$simcatprob.parm
  if(!is.null(att.profiles)){
    if(ncol(att.profiles) != K){stop("Number of attributes in Q and in att.profiles must be the same.")}
  }

  #---------------
  # Generate data
  #---------------

  simdat <- array(NA, dim = c(N, J, R))
  simalpha <- array(NA, dim = c(N, K, R))
  for(r in 1:R){
    if(!is.null(att.profiles)){
      tmp <- GDINA::simGDINA(N, Q, catprob.parm = catprob.parm, attribute = att.profiles)
      simalpha[,,r] <- att.profiles
    } else {
      switch(att.dist,
             uniform = tmp <- GDINA::simGDINA(N, Q, catprob.parm = catprob.parm, att.dist = "uniform"),
             mvnorm = tmp <- GDINA::simGDINA(N, Q, catprob.parm = catprob.parm, att.dist = "mvnorm",
                                             mvnorm.parm = list(mean = mvnorm.parm$mean,
                                                                sigma = mvnorm.parm$sigma,
                                                                cutoffs = mvnorm.parm$cutoffs)),
             higher.order = tmp <- GDINA::simGDINA(N, Q, catprob.parm = catprob.parm, att.dist = "higher.order",
                                                   higher.order.parm = list(theta = higher.order.parm$theta,
                                                                            lambda = higher.order.parm$lambda)),
             categorical = tmp <- GDINA::simGDINA(N, Q, catprob.parm = catprob.parm, att.dist = "categotical",
                                                  att.prior = categorical.parm$att.prior)
      )
      simalpha[,,r] <- tmp$attribute
    }
    simdat[,,r] <- tmp$dat
  }
  if(R == 1){
    simdat <- simdat[,,1]
    simalpha <- simalpha[,,1]
  }

  #----------------
  # Return results
  #----------------

  res <- list()
  res$simdat <- simdat
  res$simQ <- Q
  res$simalpha <- simalpha
  res$specifications <- list(N = N, R = R, item.bank = item.bank, att.profiles = att.profiles,
                             att.dist = att.dist, mvnorm.parm = mvnorm.parm, higher.order.parm = higher.order.parm,
                             categorical.parm = categorical.parm, seed = seed)
  return(res)
}
