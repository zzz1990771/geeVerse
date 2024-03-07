#' Generate Data for Simulation
#'
#' This function generates simulated data including the predictor matrix `X` and the response vector `y`,
#' based on the specified parameters. The function allows for the simulation of data under different settings
#' of correlation, distribution, and the number of observations and subjects.
#'
#' @param n_sub Integer, the number of subjects.
#' @param n_obs Integer or numeric vector, the number of observations per subject.
#' @param p Integer, the number of predictors.
#' @param beta0 Numeric vector, initial coefficients for the first few predictors.
#' @param rho Numeric, the correlation coefficient used in generating correlated errors.
#' @param type Character, the type of correlation structure (default is autoregressive).
#' @param dis Character, the distribution of errors ("normal" or "t").
#' @param ka 1 for heterogeneous errors and 0 for homogeneous errors.
#' @return A list containing two elements: `X`, the matrix of predictors, and `y`, the response vector.
#' @examples
#' sim_data <- generateData(n_sub = 100, n_obs = rep(10, 100),  p = 200,
#'                          beta0 = rep(1,7), rho = 0.6, type = "ar",
#'                           dis = "normal", ka = 1)
#' @export
generateData <- function(n_sub, n_obs, p, beta0, rho, type = "ar", dis = "normal" , ka) {
  p0 <- length(beta0)
  beta <- c(beta0, rep(0, p - p0))
  e <- NULL
  id <- NULL

  for (i in 1:n_sub) {
    id <- c(id, rep(i, n_obs[i]))
    sigmai <- Siga_cov(rho, type, n_obs[i])
    if (dis == "normal") ei <- mvtnorm::rmvnorm(1, mean = rep(0, n_obs[i]), sigma = sigmai)
    if (dis == "t") ei <- mvtnorm::rmvt(1, sigma = sigmai, df = 4, delta = rep(0, n_obs[i]))
    e <- c(e, ei)
  }

  N <- sum(n_obs)
  x <- matrix(stats::rnorm(N * p), N, p)
  y <- x %*% beta + (1 + ka * abs(x[, 1])) * e

  return(list(X = x, y = y))
}
