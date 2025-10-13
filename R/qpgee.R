#' Quantile Penalized Generalized Estimating Equations (QPGEE)
#'
#' @description Fits a quantile penalized generalized estimating equation (QPGEE)
#' model for longitudinal data using penalized quantile regression with
#' different working correlation structures.
#'
#' @param x A formula or a matrix of predictors.
#' @param ... Other arguments passed to methods.
#'
#' @return An object of class `qpgee`.
#'
#' @examples
#' # Quick Example:
#' # 1. Generate some data
#' set.seed(123)
#' sim_data <- generate_data(
#'   nsub = 50, nobs = rep(5, 50), p = 10,
#'   beta0 = c(rep(1, 5), rep(0, 5)), rho = 0.3
#' )
#'
#' # 2. Fit the model using the formula interface
#' fit <- qpgee(
#'   y ~ . - id,
#'   data = sim_data,
#'   id = sim_data$id,
#'   tau = 0.5,
#'   method = "HBIC"
#' )
#'
#' # 3. View the summary of the results
#' summary(fit)
#'
#' @export
qpgee <- function(x, ...) {
  UseMethod("qpgee")
}

#' @export
#' @rdname qpgee
#' @param id A vector identifying the clusters (subjects).
#' @param data An optional data frame.
#' @importFrom stats model.frame model.response
qpgee.formula <- function(x, id, data = parent.frame(), ...) {
  # Rename 'x' back to 'formula' for clarity inside the function
  formula <- x

  # Capture the call
  call <- match.call()

  # --- 2. Construct and Evaluate a Call to model.frame() ---
  mf_call <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "id", "data", "subset", "na.action"), names(mf_call), 0L)
  mf_call <- mf_call[c(1L, m)]
  mf_call$drop.unused.levels <- TRUE
  mf_call[[1L]] <- as.name("model.frame")

  # Cluster information
  id_var <- eval(substitute(id), envir = data, enclos = parent.frame())
  if (is.null(id_var)) {
    stop("'id' variable not found in the data frame or the environment.")
  }
  nobs <- as.vector(table(id_var))

  # --- 3. Extract Components from the Model Frame ---
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf)

  # Call the default method
  fit <- qpgee.default(x = x, y = y, nobs = nobs, ...)

  # Add call and formula to the fit object
  fit$call <- call
  fit$formula <- formula
  fit$terms <- mt
  fit$contrasts <- attr(x, "contrasts")
  fit$data <- data

  return(fit)
}


#' @param x A matrix of predictors.
#' @param y A numeric vector of response variables.
#' @param nobs A numeric vector of observations per subject.
#' @param tau The quantile to be estimated (default is 0.5).
#' @param corstr A string specifying the working correlation structure.
#'        Options include "exchangeable" (Exchangeable), "AR1" (Autoregressive),
#'        "Tri" (Tri-diagonal), "independence" (Independent), and "unstructured".
#' @param lambda A vector of penalty parameters. If NULL, auto-selection is performed.
#' @param method Criterion for penalty selection ("HBIC" or "CV").
#' @param intercept Logical; if TRUE, an intercept is added.
#' @param betaint Initial values for the beta coefficients. If NULL,
#'        non-longitudinal quantile regression is used for initialization.
#' @param nfold The number of folds used in cross-validation.
#' @param ncore Number of cores for parallel processing.
#' @param control A list of control parameters from `qpgeeControl()`, such as
#'        max_it, epsilon, shrinkCutoff, standardize and trace.
#'
#' @export
#' @rdname qpgee
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar% %do% foreach
#' @importFrom parallel stopCluster makeCluster detectCores
#' @export
qpgee.default <- function(x, y, nobs, tau = 0.5, corstr = "exchangeable",
                          lambda = NULL, method = "HBIC", intercept = TRUE,
                          betaint = NULL, nfold = 5,
                          ncore = 1, control = qpgeeControl(), ...) {

  # --- Input Validation ---
  # Validate and normalize the corstr argument to ensure it's one of the allowed values.
  corstr <- match.arg(corstr, c("exchangeable", "AR1", "Tri", "independence", "unstructured"))
  intercept = attr(x, "assign")[1] == 0
  if (!is.matrix(x)) x <- as.matrix(x)
  n_vars <- ncol(x)
  n_subjects <- length(nobs)

  args <- list(...)
  penalty_factor <- args$penalty_factor
  if (is.null(penalty_factor)) {
    penalty_factor <- rep(1, n_vars)
  }

  if (is.null(lambda)) {
    lambda_max <- 10
    lambda_min_ratio <- ifelse(n_subjects > n_vars, 1e-3, 1e-2)
    lambda <- exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio), length.out = 30))
  }

  if (length(lambda) > 1) {
    is_parallel <- (ncore > 1)
    if (is_parallel) {
      cl <- parallel::makeCluster(min(ncore, parallel::detectCores()))
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl))
    }
    `%loop%` <- if (is_parallel) `%dopar%` else `%do%`

    if (method == "HBIC") {
      fit_results <- foreach(
        l = lambda, .combine = rbind,
        .packages = c("MASS", "geeVerse")
      ) %loop% {
        res <- tryCatch({
          fit <- .fit_single_lambda(x, y, nobs, tau, corstr, l, penalty_factor,
                                    intercept, betaint, control, ...)
          list(hbic = fit$hbic)
        }, error = function(e) { list(hbic = Inf) })
        c(lambda = l, hbic = res$hbic)
      }

      fit_results <- as.data.frame(fit_results)
      best_lambda <- fit_results$lambda[which.min(fit_results$hbic)]

    } else if (method == "CV") {
      l <- NULL
      fit_results <- foreach(
        l = lambda, .combine = rbind,
        .packages = c("MASS", "geeVerse")
      ) %loop% {
        res <- tryCatch({
          subject_id_vector <- rep(1:n_subjects, times = nobs)
          nfold <- 5
          fold_assignments <- sample(rep(1:nfold, length.out = n_subjects))
          fold_errors <- numeric(nfold)
          for (k in 1:nfold) {
            test_subjects <- which(fold_assignments == k)
            train_subjects <- which(fold_assignments != k)

            test_indices <- which(subject_id_vector %in% test_subjects)
            train_indices <- which(subject_id_vector %in% train_subjects)

            # Fit model on training data
            train_fit <- .fit_single_lambda(
              x = x[train_indices, , drop = FALSE],
              y = y[train_indices],
              nobs = nobs[train_subjects],
              tau = tau, corstr = corstr, lambda = l,
              penalty_factor = penalty_factor, intercept = intercept,
              betaint = betaint, control = control, ...
            )

            # Predict on test data
            x_test <- x[test_indices, , drop = FALSE]
            preds <- as.vector(x_test %*% train_fit$coefficients)

            # Calculate check loss on test data
            fold_errors[k] <- mean(check_loss(y[test_indices] - preds, tau))
          }
          list(cv_loss = mean(fold_errors))
        }, error = function(e) { list(cv_loss = Inf) })

        c(lambda = l, cv_loss = res$cv_loss)
      }
      fit_results <- as.data.frame(fit_results)
      best_lambda <- fit_results$lambda[which.min(fit_results$cv_loss)]

    } else {
      stop("Method must be 'HBIC' or 'CV'")
    }

    if (length(best_lambda) > 1) best_lambda <- best_lambda[1]

  } else {
    best_lambda <- lambda[1]
  }

  final_fit <- .fit_single_lambda(x, y, nobs, tau, corstr, best_lambda,
                                  penalty_factor, intercept, betaint, control, ...)

  final_fit$lambda_path <- if(exists("fit_results")) fit_results else NULL
  final_fit$best_lambda <- best_lambda
  final_fit$call <- match.call()
  final_fit$n_obs <- sum(nobs)
  final_fit$n_subjects <- n_subjects

  class(final_fit) <- "qpgee"
  return(final_fit)
}


# This function is not exported. It's a helper for qpgee.default.
# It contains the core fitting algorithm from your original qpgee.est
#' @importFrom stats na.omit sd
.fit_single_lambda <- function(x, y, nobs, tau, corstr, lambda, penalty_factor,
                              intercept, betaint, control) {
  # --- 1. Initialization & Preprocessing ---
  # Extract control parameters
  max_it <- control$maxit
  epsilon <- control$epsilon
  cutoff <- control$shrinkCutoff

  # Handle intercept based on `intercept`

  if (intercept) {
    x_all <- x
    # Adjust penalty factor for intercept
    current_penalty_factor <- c(0, penalty_factor)
  } else {
    x_all <- x
    current_penalty_factor <- penalty_factor
  }

  # --- Standardization Logic ---
  was_standardized <- FALSE
  if (lambda > 0 && !is.null(control$standardize) && control$standardize) {
    # Identify columns to exclude from standardization (e.g., binary, SNP)
    cols_to_exclude <- apply(x_all, 2, function(col) {
      unique_vals <- unique(stats::na.omit(col))
      # Exclude intercept, binary, or SNP-like (0,1,2) variables
      (length(unique_vals) == 1 && unique_vals[1] == 1) || # Intercept
        (length(unique_vals) <= 2) || # Binary
        (length(unique_vals) <= 3 && all(unique_vals %in% c(0, 1, 2))) # SNP
    })

    # Get means and sds for columns that WILL be standardized
    cols_to_standardize <- !cols_to_exclude
    X_mean <- apply(x_all[, cols_to_standardize, drop = FALSE], 2, mean)
    X_sd <- apply(x_all[, cols_to_standardize, drop = FALSE], 2, stats::sd)

    # Avoid division by zero for constant columns that were not excluded
    X_sd[X_sd < .Machine$double.eps] <- 1

    # Standardize the appropriate columns
    x_all[, cols_to_standardize] <- scale(x_all[, cols_to_standardize, drop = FALSE], center = X_mean, scale = X_sd)
    was_standardized <- TRUE

    message("Note: Predictors were standardized for model fitting.")
  }

  n_subjects <- length(nobs)
  n_vars_total <- ncol(x_all)
  cluster_indices <- c(0, cumsum(nobs))

  # Initial beta values
  if (!is.null(betaint)) {
    betaint <- betaint
  } else {
    betaint <- as.numeric(quantreg::rq.fit.br(x_all, y, tau = tau)$coefficients)
  }
  beta_current <- betaint
  x_active <- x_all
  active_set <- 1:n_vars_total # Indices of active predictors

  # Iteration loop
  iter <- 0
  step_size <- 1
  max_coef_change <- Inf

  # --- 2. Iterative Fitting Loop ---
  while (max_coef_change > epsilon && iter < max_it) {
    iter <- iter + 1

    mu <- x_active %*% beta_current

    # --- 2a. Update Working Correlation Matrix R ---
    R <- .update_correlation(y - mu, nobs, corstr, tau, cluster_indices)

    # --- 2b. Update Beta Coefficients ---
    update_result <- .update_beta(x_active, y, mu, beta_current, R, nobs, tau, n_subjects, cluster_indices, lambda, current_penalty_factor[active_set])

    beta_update <- update_result$beta_update

    # Simple line search
    # (A more advanced one could be implemented here)
    mcl_old <- mean(check_loss(y - mu, tau))
    beta_new <- beta_current + step_size * beta_update
    mcl_new <- mean(check_loss(y - (x_active %*% beta_new), tau))

    if (mcl_new > mcl_old * 1.1) {
      step_size <- step_size * control$decay
    }

    beta_current <- beta_current + step_size * beta_update

    # --- 2c. Shrink Coefficients and Update Active Set ---
    if (lambda > 0) {
      shrunk_indices <- which(abs(beta_current) < cutoff)

      if (length(shrunk_indices) > 0) {
        # Ensure protected variables are not removed
        protected <- which(current_penalty_factor[active_set] == 0)
        shrunk_indices <- setdiff(shrunk_indices, protected)

        if (length(shrunk_indices) > 0) {
          # Remove from active set
          beta_current <- beta_current[-shrunk_indices]
          x_active <- x_active[, -shrunk_indices, drop = FALSE]
          active_set <- active_set[-shrunk_indices]
        }
      }
    }

    if (ncol(x_active) == 0) break

    max_coef_change <- max(abs(step_size * beta_update))
  }

  # --- 3. Finalization & Output ---
  converged <- (max_coef_change <= epsilon && iter < max_it)

  # Reconstruct full beta vector
  beta_final <- rep(0, n_vars_total)
  names(beta_final) <- colnames(x_all)
  if(length(active_set) > 0) beta_final[active_set] <- beta_current

  # --- Back-transform Coefficients ---
  if (was_standardized) {
    # Get the names of the columns that were actually standardized
    standardized_names <- names(X_mean)

    # Back-transform the coefficients for the standardized variables
    beta_final[standardized_names] <- beta_final[standardized_names] / X_sd

    # Adjust the intercept if it exists
    if (intercept) {
      intercept_adjustment <- sum(beta_final[standardized_names] * X_mean)
      beta_final["(Intercept)"] <- beta_final["(Intercept)"] - intercept_adjustment
    }
    message("Note: Coefficients were back-transformed to the original scale.")
  }

  # Calculate final metrics
  fitted_values <- x_all %*% beta_final
  residuals <- y - fitted_values
  final_mcl <- mean(check_loss(residuals, tau))

  # HBIC calculation
  n_selected <- length(active_set)
  hbic <- log(final_mcl * n_subjects) + (log(n_subjects) / (2 * n_subjects)) * log(log(n_vars_total)) * n_selected

  # Return a structured list
  list(
    coefficients = beta_final,
    fitted.values = fitted_values,
    residuals = residuals,
    R = R,
    X_selected = active_set,
    mcl = final_mcl,
    hbic = hbic,
    converged = converged,
    was_standardized = was_standardized,
    iterations = iter,
    corstr = corstr
  )
}

#' Update the Working Correlation Matrix
#'
#' This is an internal helper function to compute the working correlation matrix
#' based on standardized residuals.
#'
#' @param residuals A numeric vector of residuals (y - mu).
#' @param nobs A numeric vector of observations per subject.
#' @param corstr A string specifying the correlation structure.
#' @param tau The quantile being estimated.
#' @param cluster_indices A numeric vector of cumulative subject observation counts.
#' @return The updated working correlation matrix `R`.
#' @keywords internal
.update_correlation <- function(residuals, nobs, corstr, tau, cluster_indices) {
  # Standardize residuals based on the sign for quantile regression
  resid_sign <- 0.5 - (residuals < 1e-10)
  sd_resid <- resid_sign / sqrt(tau * (1 - tau))
  sd_resid <- scale(sd_resid) # Standardize to have mean 0, sd 1
  n_subjects <- length(nobs)
  max_nobs <- max(nobs)
  R <- matrix(0, max_nobs, max_nobs)

  if (corstr == "independence") {
    return(diag(max_nobs))
  }

  # Calculate correlation parameter(s)
  if (corstr == "exchangeable") {
    sum_prod <- 0
    n_pairs <- 0
    for (i in 1:n_subjects) {
      if (nobs[i] < 2) next
      subject_resids <- sd_resid[(cluster_indices[i] + 1):cluster_indices[i + 1]]
      sum_prod <- sum_prod + sum(subject_resids %o% subject_resids) - sum(subject_resids^2)
      n_pairs <- n_pairs + nobs[i] * (nobs[i] - 1)
    }
    alpha <- if (n_pairs > 0) sum_prod / n_pairs else 0
    R <- matrix(alpha, max_nobs, max_nobs)
    diag(R) <- 1

  } else if (corstr == "AR1") {
    sum_lag1 <- 0
    n_lag1 <- 0
    for (i in 1:n_subjects) {
      if (nobs[i] < 2) next
      subject_resids <- sd_resid[(cluster_indices[i] + 1):cluster_indices[i + 1]]
      sum_lag1 <- sum_lag1 + sum(subject_resids[-1] * subject_resids[-nobs[i]])
      n_lag1 <- n_lag1 + nobs[i] - 1
    }
    alpha <- if (n_lag1 > 0) sum_lag1 / n_lag1 else 0
    R <- alpha^abs(row(R) - col(R))

  } else if (corstr == "Tri") {
    # Assuming Tri-diagonal is similar to AR1 but only for lag 1
    sum_lag1 <- 0
    n_lag1 <- 0
    for (i in 1:n_subjects) {
      if (nobs[i] < 2) next
      subject_resids <- sd_resid[(cluster_indices[i] + 1):cluster_indices[i + 1]]
      sum_lag1 <- sum_lag1 + sum(subject_resids[-1] * subject_resids[-nobs[i]])
      n_lag1 <- n_lag1 + nobs[i] - 1
    }
    alpha <- if (n_lag1 > 0) sum_lag1 / n_lag1 else 0
    R[abs(row(R) - col(R)) == 1] <- alpha
    diag(R) <- 1

  } else if (corstr == "unstructured") {
    sum_jk <- matrix(0, max_nobs, max_nobs)
    for (i in 1:n_subjects) {
      ni <- nobs[i]
      if (ni == 0) next
      subject_resids <- sd_resid[(cluster_indices[i] + 1):cluster_indices[i + 1]]
      R[1:ni, 1:ni] <- R[1:ni, 1:ni] + subject_resids %o% subject_resids
      sum_jk[1:ni, 1:ni] <- sum_jk[1:ni, 1:ni] + 1
    }
    # Avoid division by zero
    sum_jk[sum_jk == 0] <- 1
    R <- R / sum_jk
    # Standardize to a correlation matrix
    R_sd <- sqrt(diag(R))
    R <- t(t(R / R_sd) / R_sd) # sweep(R, 2, R_sd, "/")
    diag(R) <- 1
  }

  return(R)
}

#' Update Beta Coefficients for One Iteration
#'
#' This is an internal helper function to compute the update step for beta.
#'
#' @importFrom stats pnorm dnorm
#' @keywords internal
.update_beta <- function(x_active,
                         y,
                         mu,
                         beta_current,
                         R,
                         nobs,
                         tau,
                         n_subjects,
                         cluster_indices,
                         lambda,
                         penalty_factor_active,
                         f0 = NULL) {
  if (is.null(f0)) {
    f0 <- rep(1, length(y))
  }

  n_total_obs <- length(y)
  n_active_vars <- ncol(x_active)

  score_vector <- matrix(0, n_active_vars, 1)
  hessian_approx <- matrix(0, n_active_vars, n_active_vars)

  # Calculate smoothing components
  r2 <- sqrt(diag(x_active %*% t(x_active)))
  r2[r2 < .Machine$double.eps] <- 1e-6 # Replace zeros with a small number
  hh <- y - mu
  fr <- (pnorm(sqrt(n_subjects) * hh / r2) - (1 - tau))

  # Loop over subjects to calculate score and hessian
  for (i in 1:n_subjects) {
    idx <- (cluster_indices[i] + 1):cluster_indices[i + 1]
    ni <- nobs[i]
    if (ni == 0)
      next

    # Subset for the current subject
    Ra <- R[1:ni, 1:ni, drop = FALSE]
    xa <- x_active[idx, , drop = FALSE]
    ha <- hh[idx]
    r2a <- r2[idx]

    # Use geninv for robustness
    Ra_inv <- geninv(Ra)

    Ga <- if (ni == 1)
      f0[idx]
    else
      diag(f0[idx])
    Ha2 <- if (ni == 1)
      dnorm(sqrt(n_subjects) * ha / r2a) / r2a
    else
      diag(dnorm(sqrt(n_subjects) * ha / r2a) / r2a)

    # Add subject's contribution
    score_vector <- score_vector + t(xa) %*% Ga %*% Ra_inv %*% fr[idx]
    hessian_approx <- hessian_approx + sqrt(n_subjects) * t(xa) %*% Ga %*% Ra_inv %*% Ha2 %*% xa
  }

  # Calculate penalty matrix
  if (lambda > 0) {
    eps <- 1e-6 # To avoid division by zero
    # SCAD penalty derivative
    scad_deriv <- pp_scad_lin(abs(as.vector(beta_current)), lambda) / (abs(as.vector(beta_current)) + eps)

    # Apply penalty factor
    pen_weights <- scad_deriv * penalty_factor_active

    # Penalty matrix
    sigma <- n_subjects * diag(pen_weights, nrow = n_active_vars)

    # Calculate update with penalty
    hessian_inv <- geninv(hessian_approx + sigma)
    beta_update <- hessian_inv %*% (score_vector - sigma %*% beta_current)
  } else {
    # Calculate update without penalty
    hessian_inv <- geninv(hessian_approx)
    beta_update <- hessian_inv %*% score_vector
  }
  return(list(beta_update = beta_update))
}

#' Control Parameters for qpgee
#'
#' @description
#' Provides control parameters for the Quantile Penalized Generalized Estimating
#' Equations (QPGEE) fitting procedure. Similar in spirit to `geese.control()` in
#' the geepack package.
#'
#' @param epsilon Convergence tolerance for the parameter estimates. Iteration stops
#'   when the maximum change in coefficients is below this value.
#' @param maxit Maximum number of iterations.
#' @param decay Decay rate of learning step.
#' @param trace Logical indicating if output should be produced for each iteration.
#'   (You can decide how much information to show inside the C++/R loop.)
#' @param standardize Logical indicating whether to scale X.
#' @param shrinkCutoff Threshold below which coefficients are shrunk to zero
#'   (removal of “small” coefficients).
#'
#' @return A list with the components \code{epsilon}, \code{maxit},
#'   \code{trace}, and \code{shrinkCutoff}.
#'
#' @examples
#' ctrl <- qpgeeControl(epsilon = 1e-5, maxit = 200, trace = TRUE)
#' @export
qpgeeControl <- function(epsilon = 1e-4,
                         decay  = 1,
                         maxit = 100,
                         trace = FALSE,
                         standardize = FALSE,
                         shrinkCutoff = 1e-4) {
  if (!is.numeric(epsilon) || epsilon <= 0)
    stop("value of 'epsilon' must be > 0")
  if (!is.numeric(maxit) || maxit <= 0)
    stop("maximum number of iterations must be > 0")

  list(
    epsilon      = epsilon,
    maxit        = as.integer(maxit),
    decay        = as.numeric(decay),
    trace        = as.logical(trace),
    standardize  = as.logical(standardize),
    shrinkCutoff = shrinkCutoff
  )
}
