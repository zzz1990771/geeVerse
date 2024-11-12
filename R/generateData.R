#' Generate Data for Simulation
#'
#' This function generates simulated data including the predictor matrix `X` and the response vector `y`,
#' based on the specified parameters. The function allows for the simulation of data under different settings
#' of correlation, distribution, and the number of observations and subjects.
#'
#' @param nsub Integer, the number of subjects.
#' @param nobs Integer or numeric vector, the number of observations per subject.
#' @param p Integer, the number of predictors.
#' @param beta0 Numeric vector, initial coefficients for the first few predictors.
#' @param rho Numeric, the correlation coefficient used in generating correlated errors.
#' @param correlation Character, the correlation of correlation structure (default is autoregressive).
#' @param dis Character, the distribution of errors ("normal" or "t").
#' @param ka 1 for heterogeneous errors and 0 for homogeneous errors.
#' @param SNPs User can provide simulated or real SNPs for genetic data simulation.
#' @return A list containing two elements: `X`, the matrix of predictors, and `y`, the response vector.
#' @examples
#' sim_data <- generateData(
#'   nsub = 100, nobs = rep(10, 100), p = 200,
#'   beta0 = c(rep(1, 7), rep(0, 193)), rho = 0.6, correlation = "AR1",
#'   dis = "normal", ka = 1
#' )
#' @export
generateData <- function(nsub, nobs, p, beta0, rho, correlation = "AR1", dis = "normal", ka = 0, SNPs = NULL) {
  beta <- beta0
  e <- NULL
  id <- NULL

  for (i in 1:nsub) {
    id <- c(id, rep(i, nobs[i]))
    sigmai <- Siga_cov(rho, correlation, nobs[i])
    if (dis == "normal") ei <- mvtnorm::rmvnorm(1, mean = rep(0, nobs[i]), sigma = sigmai)
    # print(length(ei))
    if (dis == "t") ei <- mvtnorm::rmvt(1, sigma = sigmai, df = 4, delta = rep(0, nobs[i]))
    e <- c(e, ei)
  }
  # print(length(e))
  N <- sum(nobs)
  if (!is.null(SNPs)) {
    SNPslon <- SNPs[rep(1:nrow(SNPs), each = 5), ]
    x <- matrix(stats::rnorm(N * p / 2), N, p / 2)
    x <- cbind(x[, 1:5], SNPslon, (x[, 6:(p / 2)]))
    colnames(x)[c(1:5, (p / 2 + 6):p)] <- paste0("control", 1:(p / 2))
  } else {
    x <- matrix(stats::rnorm(N * p), N, p)
  }

  y <- as.vector(x %*% beta + (1 + ka * abs(x[, 1])) * e)

  return(list(X = x, y = y))
}
