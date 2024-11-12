#' GeeVerse: Wrapper for Quantile Penalized Generalized Estimating Equations
#'
#' This function is a wrapper for qpgee that allows running the model for multiple quantile levels (tau).
#'
#' @param x A matrix of predictors.
#' @param y A numeric vector of response variables.
#' @param tau A vector of quantiles to be estimated (default is c(0.25, 0.5, 0.75)).
#' @param ... Additional arguments to be passed to qpgee function.
#' @return A list containing the results for each tau value and a combined beta matrix.
#' @export
#' @importFrom foreach %do% foreach
geeVerse <- function(x, y, tau = c(0.25, 0.5, 0.75), ...) {
  # Ensure tau is a vector
  if (!is.vector(tau)) {
    stop("tau must be a vector of quantile values")
  }

  # Run qpgee for each tau
  results <- foreach::foreach(t = tau) %do% {
    qpgee(x = x, y = y, tau = t, ...)
  }

  # Name the results list with tau values
  names(results) <- paste0("tau_", tau)

  # Compile beta estimates across all tau values
  beta_matrix <- do.call(cbind, lapply(results, function(r) r$beta))
  colnames(beta_matrix) <- paste0("tau_", tau)

  # Add the combined beta matrix to the results
  results$beta_matrix <- beta_matrix

  return(results)
}
