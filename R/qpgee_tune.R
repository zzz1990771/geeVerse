#' Quantile Penalized Generalized Estimating Equations with Auto Selected Penalty level
#'
#' This function automatically select the penalty level by going through a list
#' of lambdas, and select the best level of penalty with high-dimensional BIC
#' (HBIC) or cross-validation (CV).
#'
#' @param x A matrix of predictors.
#' @param y A numeric vector of response variables.
#' @param tau The quantile to be estimated (default is 0.5, the median).
#' @param method The criterion to select level of penalty. Currently it only
#'        supports "HBIC".
#' @param ncore A numeric value specifying how many core to use.
#' @param nk A numeric vector indicating the number of observations per subject.
#' @param worktype A string specifying the working correlation structure.
#'        Options include "CS" (Compound Symmetry), "AR" (Autoregressive),
#'        "Tri" (Tri-diagonal), and "Ind" (Independent).
#' @param lambda A vector of penalty parameter for regularization. If not provided,
#' a grid will be provided by this function.
#' @param betaint Initial values for the beta coefficients. If NULL,
#'        non-longitudinal quantile regression is used for initialization.
#' @param f0 estimated conditional error distributions.
#' @param max_it Maximum number of iterations (default is 100).
#' @param cutoff Threshold for coefficient shrinkage (default is 0.1).
#' @return A list containing the following components:
#'           \item{beta}{Estimated beta coefficients.}
#'           \item{g}{Fitted values of the linear predictor.}
#'           \item{R}{Estimated working correlation matrix.}
#'           \item{X_selected}{Indices of selected predictors.}
#'           \item{mcl}{Mean check loss.}
#'           \item{hbic}{Hannan-Quinn Information Criterion value.}
#'           \item{converge}{Boolean indicating whether the algorithm converged.}
#' @examples
#' # Example usage:
#'
#' sim_data <- generateData(n_sub = 20, n_obs = rep(10, 20),  p = 20,
#'                          beta0 = rep(1,5), rho = 0.1, type = "ar",
#'                           dis = "normal", ka = 1)
#'
#' X=sim_data$X
#' y=sim_data$y
#'
#' #fit qpgee with auto selected lambda
#' qpgee.fit = qpgee_tune(X,y,tau=0.5,nk=rep(10, 20),ncore=1)
#' qpgee.fit$beta
#'
#' @export
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar% %do% foreach
#' @importFrom parallel stopCluster makeCluster detectCores
qpgee_tune <-
  function(x, y,
           tau = 0.5,
           method = "HBIC",
           ncore = 1,
           nk = rep(1, length(y)),
           worktype = "CS",
           lambda = NULL,
           f0 = NULL,
           betaint = NULL,
           max_it = 100,
           cutoff = 10 ^ -1) {
    #setting up lambda
    if (is.null(lambda)) {
      #similiar to glmnet
      lambda_max = 10
      lambda.min.ratio = ifelse(length(nk) > NCOL(x), 1e-03, 0.01)
      lambda = exp(seq(
        log(lambda_max),
        log(lambda_max * lambda.min.ratio),
        length.out = 30
      ))
    }

    if (method == "HBIC"){
      if (ncore == 1) {
        l = lambda
        fit_results <-foreach::foreach(l = l, .combine = rbind,
                                       .packages = c("MASS","geeVerse")) %do% {
                                         result <-qpgee(x, y, tau,
                                                        nk, worktype,
                                                        l, betaint, f0,
                                                        max_it, cutoff)
                                         c(l, result$mcl, result$hbic,
                                           paste0(result$X_selected, collapse = ""))
                                       }
      }else if (ncore > 1 && ncore%%1 == 0){
        cl <- makeCluster(min(ncore, detectCores()))
        registerDoParallel(cl)
        fit_results <-foreach::foreach(l = lambda, .combine = rbind,
                                       .packages = c("MASS","geeVerse")) %dopar% {
                                         result <-qpgee(x, y, tau,
                                                        nk, worktype,
                                                        l, betaint, f0,
                                                        max_it, cutoff)
                                         c(l, result$mcl, result$hbic,
                                           paste0(result$X_selected, collapse = ""))
                                       }
        parallel::stopCluster(cl)
      }else{
        stop("ncore must be a integer greater than 1")
      }
      fit_results = as.data.frame(fit_results)
      fit_results[, 3] = as.numeric(fit_results[, 3])
      best_lambda <-
        lambda[which(fit_results[, 3] == min(fit_results[, 3]))]
      if (length(best_lambda) > 1) {
        best_lambda = best_lambda[1]
      }
    }else{
      stop("CV is currently in development")
    }
    fit = qpgee(x, y, tau,
          nk, worktype,
          best_lambda, betaint, f0,
          max_it, cutoff)
    fit$best_lambda = best_lambda
    return(fit)
  }

# qpgee_cv <-
#   function(x,
#            y,
#            betaint,
#            nk,
#            worktype,
#            f0,
#            tau = 0.5,
#            lambda = NULL,
#            max_it = 100,
#            cutoff = 10 ^ -3,
#            nfold = 3) {
#     #setting up lambda
#     if (is.null(lambda)) {
#       #similiar to glmnet
#       lambda_max = 10
#       lambda.min.ratio = ifelse(length(nk) > NCOL(x), 1e-03, 0.01)
#       lambda = exp(seq(
#         log(lambda_max),
#         log(lambda_max * lambda.min.ratio),
#         length.out = 30
#       ))
#     }
#     for (l in lambda) {
#       #this method only applies to balanced data, nk are the same
#       n_obs = nk[1]
#       n_sub = length(nk)
#       sampled_sind = suppressWarnings(split(sample(1:n_sub), 1:nfold))
#       cl = c()
#       for (fold in 1:nfold) {
#         test_sind = sampled_sind[[fold]]
#         train_sind = setdiff(1:n_sub, test_sind)
#         test_ind = as.vector(sapply((test_sind - 1) * n_obs + 1, function(x)
#           x:(x + n_obs - 1)))
#         train_ind = as.vector(sapply((train_sind - 1) * n_obs + 1, function(x)
#           x:(x + n_obs - 1)))
#
#         result <-
#           qlingee_scad(x[train_ind, ],
#                        y[train_ind],
#                        betaint,
#                        rep(n_obs, length(train_sind)),
#                        worktype,
#                        f0,
#                        tau,
#                        l,
#                        max_it,
#                        cutoff)
#         #use check loss to select best lambda via cross-validation
#         cl = c(cl, check_loss(y[test_ind] - x[test_ind, ] %*% result$beta, tau))
#       }
#
#       mcl = mean(cl)
#       print(mcl)
#     }
#
#   }



