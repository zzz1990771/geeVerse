# ---
# Example Script to Test the Refactored `qpgee` Package
# ---

# Step 1: Load Required Libraries
# ---------------------------------
# Ensure you have the necessary packages installed.
# install.packages("mvtnorm")
library(mvtnorm)

# Load your package functions here. For example:
library(geeVerse)
# Or source the files directly if they are not in a package yet:
# source("path/to/qpgee_functions.R")


# Step 2: Define Helper and Data Generation Functions
# ---------------------------------------------------

#' Create a Covariance Matrix (Helper for generate_data)
#'
#' This function was a missing dependency for `generate_data`. It creates a
#' covariance matrix based on a given correlation structure.
#'
#' @param rho The correlation coefficient.
#' @param corstr The correlation structure ("AR1", "exchangeable").
#' @param n The dimension of the matrix.
#' @return A covariance matrix.
Siga_cov <- function(rho, corstr, n) {
  if (corstr == "AR1") {
    # Autoregressive structure
    m <- abs(outer(1:n, 1:n, "-"))
    return(rho^m)
  } else if (corstr == "exchangeable") {
    # Exchangeable (compound symmetric) structure
    m <- matrix(rho, n, n)
    diag(m) <- 1
    return(m)
  } else {
    # Default to independence
    return(diag(n))
  }
}


#' Generate Data for Simulation
#'
#' This function generates simulated data including the predictor matrix `X` and the response vector `y`,
#' based on the specified parameters.
#'
#' @param nsub Integer, the number of subjects.
#' @param nobs Integer or numeric vector, the number of observations per subject.
#' @param p Integer, the number of predictors.
#' @param beta0 Numeric vector, the true coefficients.
#' @param rho Numeric, the correlation coefficient for errors.
#' @param corstr Character, the correlation structure ("AR1", "exchangeable").
#' @param dis Character, the distribution of errors ("normal" or "t").
#' @param ka Heterogeneity parameter for errors.
#' @return A data frame containing `y`, `id`, and predictors.
generate_data <- function(nsub, nobs, p, beta0, rho, corstr = "AR1", dis = "normal", ka = 0) {
  beta <- beta0
  e <- NULL
  id <- NULL

  for (i in 1:nsub) {
    id <- c(id, rep(i, nobs[i]))
    sigmai <- Siga_cov(rho, corstr, nobs[i])
    if (dis == "normal") ei <- mvtnorm::rmvnorm(1, mean = rep(0, nobs[i]), sigma = sigmai)
    if (dis == "t") ei <- mvtnorm::rmvt(1, sigma = sigmai, df = 4, delta = rep(0, nobs[i]))
    e <- c(e, ei)
  }

  N <- sum(nobs)
  x <- matrix(stats::rnorm(N * p), N, p)
  colnames(x) <- paste0("x", 1:p)

  y <- as.vector(x %*% beta + (1 + ka * abs(x[, 1])) * e)
  sim_data <- as.data.frame(cbind(y, id, x))

  return(sim_data)
}


# Step 3: Set Simulation Parameters and Generate Data
# ---------------------------------------------------
cat("--- Generating Simulation Data ---\n")

# Parameters for a moderately-sized problem
n_subjects <- 50
obs_per_subject <- rep(5, n_subjects) # 5 observations per subject
n_predictors <- 40

# Define true coefficients: 5 non-zero, the rest are zero
true_beta <- c(rep(1.5, 3), rep(-1.5, 2), rep(0, n_predictors - 5))

# Generate the data
set.seed(123) # for reproducibility
sim_data <- generate_data(
  nsub = n_subjects,
  nobs = obs_per_subject,
  p = n_predictors,
  beta0 = true_beta,
  rho = 0.5,
  corstr = "AR1",
  dis = "normal"
)

cat("Data generated with", n_subjects, "subjects and", n_predictors, "predictors.\n\n")


# Step 4: Fit the QPGEE Model using the Formula Interface
# -------------------------------------------------------
cat("--- Fitting QPGEE Model ---\n")
cat("This may take a moment, especially with lambda selection...\n")

# Use the formula interface, which is the recommended way for users.
# `y ~ . - id` means use y as the response and all other columns as predictors, except for 'id'.
# We set ncore = 2 to demonstrate parallel processing for lambda selection.
# If you have more cores, you can increase this number.
fit <- qpgee(
  formula = y ~ . - id,
  data = sim_data,
  id = sim_data$id,
  tau = 0.1, # Median regression
  corstr = "exchangeable",
  method = "HBIC", # Use HBIC to select the best lambda
  ncore = 10, # Use 2 cores for parallel computation
  intercept = FALSE,
  lambda = NULL,
  control = qpgeeControl(maxit = 100, standardize = TRUE)
)

fit_old <- qpgee_old(
  formula = y ~ . - id,
  data = sim_data,
  id = sim_data$id,
  tau = 0.1, # Median regression
  corstr = "exchangeable",
  method = "HBIC", # Use HBIC to select the best lambda
  ncore = 10, # Use 2 cores for parallel computation
  intercept = FALSE,
  lambda = NULL,
  control = qpgeeControl(maxit = 100)
)
cbind(fit$coefficients,fit_old$beta)

cat("Model fitting complete.\n\n")


# Step 5: Examine the Results using S3 Methods
# --------------------------------------------

# --- 5a. Basic Print Output ---
cat("--- 1. Basic Print Output (`print(fit)`) ---\n")
print(fit)
cat("\n\n")

# --- 5b. Detailed Summary ---
cat("--- 2. Detailed Summary (`summary(fit)`) ---\n")
model_summary <- summary(fit)
print(model_summary)
cat("\n\n")

# --- 5c. Extract Coefficients ---
cat("--- 3. Extracting Coefficients (`coef(fit)`) ---\n")
coefficients <- coef(fit)
# Show only the non-zero coefficients for clarity
non_zero_coefs <- coefficients[coefficients != 0]
print(non_zero_coefs)
cat("\nComparison with true non-zero coefficients:\n")
print(true_beta[1:5])
cat("\n\n")


# Step 6: Make Predictions on New Data
# --------------------------------------
cat("--- 4. Making Predictions (`predict(fit, newdata)`) ---\n")

# Create a small new dataset with the same structure
set.seed(456)
new_data <- generate_data(
  nsub = 2,
  nobs = rep(3, 2), # 3 observations for 2 new subjects
  p = n_predictors,
  beta0 = true_beta,
  rho = 0.5,
  corstr = "AR1"
)

# Use the predict method
predictions <- predict(fit, newdata = new_data)

cat("Generated predictions for 2 new subjects (6 total observations):\n")
print(predictions)
cat("\n--- Test Script Finished ---\n")

