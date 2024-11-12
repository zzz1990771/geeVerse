# scad penalty function
pp_scad_sim <- function(z, lambda_si, lambda_pl = 0, a = 3.7) {
  # sim model so dsi always =length(z)
  dsi <- length(z)
  x <- matrix(abs(z), ncol = 1)
  value <- matrix(c(rep(0, dim(x)[[1]])), ncol = 1)
  # penalize separately here
  for (i in 1:dim(x)[[1]]) {
    lambda <- ifelse(i <= dsi, lambda_si, lambda_pl)
    # print(sprintf("coeff is %g and lambda is %g",i,lambda))

    if (x[i] < lambda) {
      value[i] <- lambda
    } else if (x[i] < a * lambda) {
      value[i] <- (a * lambda - x[i]) / (a - 1)
    } else {
      value[i] <- 0
    }
  }
  return(value)
}


# correlated error function
#' Generate Covariance Matrix
#'
#' This function generates a covariance matrix based on the specified correlation structure.
#' The function supports "compound symmetry" (cs) and "autoregressive" (ar) correlation structures,
#' as well as an identity matrix as the default option when neither "cs" nor "AR1" is specified.
#'
#' @param rho Numeric, the correlation coefficient used for generating the covariance matrix.
#'        For "cs" or "exchangeable", it represents the common correlation between any two observations.
#'        For "AR1", it represents the correlation between two consecutive observations,
#'        with the correlation decreasing for observations further apart.
#' @param correlation Character, specifies the correlation of correlation structure for the covariance matrix.
#'        Options are "cs" or "exchangeable" for compound symmetry, "AR1" for autoregressive, and any other input
#'        will result in an identity matrix.
#' @param nt Integer, the dimension of the square covariance matrix (number of time points or observations).
#'
#' @return A square matrix of dimension `nt` representing the specified covariance structure.
#'
#' @export
Siga_cov <- function(rho, correlation, nt) {
  sigma <- matrix(0, nt, nt)
  if (correlation == "cs" || correlation == "exchangeable") {
    sigma <- (1 - rho) * diag(nt) + rho * matrix(1, nt, nt)
  } else if (correlation == "AR1") {
    for (i in 1:nt) {
      for (j in 1:nt) {
        sigma[i, j] <- rho^(abs(i - j))
      }
    }
  } else {
    sigma <- diag(nt)
  }
  return(sigma)
}

# check loss
check_loss <- function(u, tau) {
  tol <- 1e-10
  u * (tau - (u < tol))
}

# log space between 10^a to 10^b
logspace <- function(d1, d2, n = 20) exp(log(10) * seq(d1, d2, length.out = n))


# compile results
#' Compile Results from qpgee()
#'
#' This function reports correct percentage, TP, FP, MSE and MAD from a (list of)
#' fitted qpgee model comparing to the true betas.
#'
#' @param qpgee_results A (list of) fitted qpgee model.
#' @param beta0 True beta used in true data generation process.
#' @param threshold Integer, the threshold to determine whether a esimated beta should
#' be consider as 0.
#'
#' @return a vector contains correct percentage, TP, FP, MSE and MAD and its standard
#' error if Monte Carlo simulations.
#'
#' @export
compile_result <- function(qpgee_results, beta0, threshold = 10^-3) {
  UseMethod("compile_result")
}

# compile results
#' Compile Results from qpgee()
#'
#' This function reports correct percentage, TP, FP, MSE and MAD from a (list of)
#' fitted qpgee model comparing to the true betas.
#'
#' @param qpgee_results A (list of) fitted qpgee model.
#' @param beta0 True beta used in true data generation process.
#' @param threshold Integer, the threshold to determine whether a esimated beta should
#' be consider as 0.
#'
#' @return a vector contains correct percentage, TP, FP, MSE and MAD and its standard
#' error if Monte Carlo simulations.
#'
#' @export
compile_result.qpgee <- function(qpgee_results, beta0, threshold = 10^-3) {
  p <- length(beta0)
  if (inherits(qpgee_results, "qpgee")) {
    qpgee_results <- list(qpgee_results)
  }
  n_sim <- length(qpgee_results)

  result_structure <- matrix(NA, nrow = n_sim, ncol = 6)
  for (sim in 1:n_sim) {
    # print(sim)
    selected <- rep(FALSE, p)
    selected <- abs(qpgee_results[[sim]]$beta) > threshold
    selected_true <- abs(beta0) > threshold
    correct <- FALSE
    TP <- sum(selected & selected_true)
    FP <- sum(selected & !selected_true)
    FN <- sum(!selected & selected_true)
    F1 <- (2 * TP) / (2 * TP + FP + FN)
    if ((TP - FP) == sum(abs(beta0) > 0)) {
      correct <- TRUE
    }
    beta_est <- rep(0, p)
    beta_est <- qpgee_results[[sim]]$beta
    MSE <- sum((beta_est - beta0)^2)
    MAD <- sum(abs(beta_est - beta0))
    result_structure[sim, ] <- c(F1, correct, TP, FP, MSE, MAD)
  }
  if (n_sim == 1) {
    results_table <- rbind(apply(result_structure, 2, mean))
  } else {
    results_table <- rbind(
      apply(result_structure, 2, mean),
      apply(result_structure, 2, stats::sd)
    )
  }
  colnames(results_table) <- c("F1", "correct%", "TP", "FP", "MSE", "MAD")
  return(results_table)
}


#' Compile Results from list of qpgee()
#'
#' This function reports correct percentage, TP, FP, MSE and MAD from a (list of)
#' fitted qpgee model comparing to the true betas.
#'
#' @param qpgee_results A (list of) fitted qpgee model.
#' @param beta0 True beta used in true data generation process.
#' @param threshold Integer, the threshold to determine whether a esimated beta should
#' be consider as 0.
#'
#' @return a vector contains correct percentage, TP, FP, MSE and MAD and its standard
#' error if Monte Carlo simulations.
#'
#' @export
compile_result.default <- function(qpgee_results, beta0, threshold = 0.1) {
  p <- length(beta0)
  if (inherits(qpgee_results, "qpgee")) {
    qpgee_results <- list(qpgee_results)
  }
  n_sim <- length(qpgee_results)

  result_structure <- matrix(NA, nrow = n_sim, ncol = 6)
  for (sim in 1:n_sim) {
    # print(sim)
    selected <- rep(FALSE, p)
    selected <- abs(qpgee_results[[sim]]$beta) > threshold
    selected_true <- abs(beta0) > threshold
    correct <- FALSE
    TP <- sum(selected & selected_true)
    FP <- sum(selected & !selected_true)
    FN <- sum(!selected & selected_true)
    F1 <- (2 * TP) / (2 * TP + FP + FN)
    if ((TP - FP) == sum(abs(beta0) > 0)) {
      correct <- TRUE
    }
    beta_est <- rep(0, p)
    beta_est <- qpgee_results[[sim]]$beta
    MSE <- sum((beta_est - beta0)^2)
    MAD <- sum(abs(beta_est - beta0))
    result_structure[sim, ] <- c(F1, correct, TP, FP, MSE, MAD)
  }
  if (n_sim == 1) {
    results_table <- rbind(apply(result_structure, 2, mean))
  } else {
    results_table <- rbind(
      apply(result_structure, 2, mean),
      apply(result_structure, 2, stats::sd)
    )
  }
  colnames(results_table) <- c("F1", "correct%", "TP", "FP", "MSE", "MAD")
  return(results_table)
}

# compile results
# compile_result_hd <- function(qpgee_results,SIS_results,beta0,threshold = 0.1){
#   if(class(qpgee_results) == "qpgee"){qpgee_results = list(qpgee_results);print("true")}
#   p = length(beta0)
#   n_sim = length(qpgee_results)
#   results_table <- matrix(NA,ncol=5,nrow=2)
#   colnames(results_table) <- c("percentage","TP","FP","MSE","MAD")
#   result_structure = matrix(NA,nrow=n_sim,ncol=5)
#   for(sim in 1:n_sim){
#     #SIS results
#     SIS_X <- SIS_results[[sim]]
#     selected <- rep(FALSE,p)
#     selected[SIS_X] <- abs(qpgee_results[[sim]]$beta)>threshold
#     selected_true <- abs(beta0)>threshold
#     correct <- FALSE
#     TP <- sum(selected & selected_true)
#     FP <- sum(selected & !selected_true)
#     if((TP-FP)==sum(abs(beta0)>0)){correct = TRUE}
#     beta_est <- rep(0,p)
#     beta_est[SIS_X] <- qpgee_results[[sim]]$beta
#     MSE <- sum((beta_est-beta0)^2)
#     MAD <- sum(abs(beta_est-beta0))
#     result_structure[sim,] <- c(correct,TP,FP,MSE,MAD)
#   }
#   results_table <- rbind(apply(result_structure,2,mean),
#                          apply(result_structure,2,sd))
#
#   return(results_table)
#
# }



cosf <- function(t, p) {
  (cos(t))^(p - 2)
}

zeroDelt <- function(t, eps, pp) {
  stopifnot(t <= 1)
  given <- eps * stats::integrate(cosf, 0, pi / 2, p = pp)$value
  stats::integrate(cosf, 0, asin(t), p = pp)$value - given
}

computeDelta <- function(eps, p) {
  # compute delta value for Singular Value upper bound. Ref: Hochstenbach (2013)
  1 / stats::uniroot(zeroDelt, interval = c(0, 0.5), eps = eps, pp = p)$root
}
