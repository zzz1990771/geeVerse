summary.qpgee1 <- function(object, ...) {
  # Extract coefficients
  coef <- object$coefficients

  # Extract variance-covariance matrix
  vcov <- object$vcov

  # Compute standard errors
  std.err <- sqrt(diag(vcov))

  # Compute Wald statistics
  wald <- (coef / std.err)^2

  # Compute p-values (assuming chi-squared distribution)
  p.value <- stats::pchisq(wald, df = 1, lower.tail = FALSE)

  # Create coefficients table
  coefficients <- data.frame(
    Estimate = coef,
    Std.Error = std.err,
    Wald = wald,
    `P-value` = p.value
  )

  # Extract correlation structure type
  corstr <- object$corstr

  # Extract correlation parameters, if they exist
  corpar <- object$corpar
  if (!is.null(corpar)) {
    # Assuming corpar is a named vector
    correlation <- data.frame(
      Parameter = names(corpar),
      Estimate = corpar
    )
  } else {
    correlation <- NULL
  }

  # Extract scale parameter, if applicable
  scale <- object$scale

  # Extract number of observations
  nobs <- object$nobs

  # Extract number of clusters, if applicable
  nclusters <- object$nclusters

  # Assemble summary object as a list
  summary_obj <- list(
    call = object$call,
    coefficients = coefficients,
    corstr = corstr,
    correlation = correlation,
    scale = scale,
    nobs = nobs,
    nclusters = nclusters
  )

  # Assign class to the summary object
  class(summary_obj) <- "summary.qpgee"

  return(summary_obj)
}


#' @export
summary.qpgee <- function(object, ...) {
  # This function now collects all necessary info for the new print method
  structure(
    list(
      formula = object$formula,
      n_subjects = object$n_subjects,
      n_obs = object$n_obs,
      coefficients = object$coefficients,
      converged = object$converged,
      best_lambda = object$best_lambda,
      hbic = object$hbic,
      corstr = object$corstr,
      R = object$R
    ),
    class = "summary.qpgee"
  )
}

#' Print summary method for qpgee model objects
#'
#' @param x A `qpgee` model object.
#' @param digits Default digits.
#' @param ... Additional arguments (not used).
#' @return Prints a summary of the qpgee model.
#' @export
print.summary.qpgee <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  # This print method is updated to match the user's requested format
  cat("Quantile Penalized Generalized Estimating Equations (QPGEE)\n")
  if (!is.null(x$formula)) {
    cat("Formula: ", deparse(x$formula), "\n")
  }
  cat("Number of Clusters: ", x$n_subjects, "\n")
  cat("Number of observations:", x$n_obs, "\n")

  cat("\nCoefficients:\n")
  if(length(x$coefficients) > 50){
    print(round(x$coefficients[1:50], digits))
    cat("...", length(x$coefficients) - 50, "coefficients are not shown.\n")
  } else {
    print(round(x$coefficients, digits))
  }

  cat("\nConverged: ", x$converged, "\n")
  cat("Best lambda: ", format(x$best_lambda, digits = digits), "\n")
  cat("HBIC: ", format(x$hbic, digits = digits), "\n")

  cat("\nCorrelation structure:", x$corstr, "\n")
  if (!is.null(x$R)) {
    cat("Estimated correlation matrix (or parameters):\n")
    # Only print the relevant part of the R matrix
    max_dim <- min(nrow(x$R), 5) # Print max 5x5
    print(round(x$R[1:max_dim, 1:max_dim], digits = digits))
  }
}

#' @export
#' @importFrom stats coef
print.qpgee <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(coef(x))
}

#' @export
coef.qpgee <- function(object, ...) {
  return(object$coefficients)
}


#' @export
#' @importFrom stats terms delete.response model.frame coef
predict.qpgee <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    return(object$fitted.values)
  }

  if (!is.data.frame(newdata)) {
    stop("'newdata' must be a data frame.")
  }

  # Ensure newdata has the right format
  tt <- terms(object$formula, "p", data = object$data) # terms object from original data
  Terms <- delete.response(tt)
  mf <- model.frame(Terms, newdata)
  X_new <- model.matrix(Terms, mf)

  # Handle intercept
  if ("(Intercept)" %in% names(coef(object))) {
    if (!("(Intercept)" %in% colnames(X_new))) {
      X_new <- cbind("(Intercept)" = 1, X_new)
    }
  }

  # Align columns
  beta <- coef(object)
  common_vars <- intersect(names(beta), colnames(X_new))
  X_aligned <- matrix(0, nrow = nrow(X_new), ncol = length(beta), dimnames = list(NULL, names(beta)))
  X_aligned[, common_vars] <- X_new[, common_vars]

  as.vector(X_aligned %*% beta)
}

