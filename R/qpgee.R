#' Quantile Penalized Generalized Estimating Equations (QPGEE)
#'
#' This function implements Quantile Penalized Generalized Estimating Equations
#' (QPGEE) for longitudinal data analysis. It estimates parameters using a
#' penalized quantile regression approach within a GEE framework, allowing for
#' different working correlation structures.
#'
#' @param x A matrix of predictors.
#' @param y A numeric vector of response variables.
#' @param tau The quantile to be estimated (default is 0.5, the median).
#' @param nk A numeric vector indicating the number of observations per subject.
#' @param worktype A string specifying the working correlation structure.
#'        Options include "CS" (Compound Symmetry), "AR" (Autoregressive),
#'        "Tri" (Tri-diagonal), and "Ind" (Independent).
#' @param lambda The penalty parameter for regularization (default is 0.1).
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
#' sim_data <- generateData(n_sub = 100, n_obs = rep(10, 100),  p = 100,
#'                          beta0 = rep(1,7), rho = 0.6, type = "ar",
#'                           dis = "normal", ka = 1)
#'
#' X=sim_data$X
#' y=sim_data$y
#'
#' #fit qpgee
#' qpgee.fit = qpgee(X,y,tau=0.5,nk=rep(10, 100))
#' qpgee.fit$beta
#'
#' @export
#' @importFrom quantreg rq
#' @importFrom stats pnorm
#' @importFrom stats dnorm
#' @useDynLib geeVerse
qpgee <-function(x, y,
                 tau = 0.5,
                 nk = rep(1, length(y)),
                 worktype = "CS",
                 lambda = 0.1,
                 betaint = NULL,
                 f0 = NULL,
                 max_it = 100,
                 cutoff = 10 ^ -1) {
  x_all = x
  cn = c(0, cumsum(nk))
  nsub = length(nk)
  nx = dim(x)[2]
  N = sum(nk)

  #if initial beta is not provided, use non-longitudinal quantile regression as initial
  if (is.null(betaint)) {
    betaint = stats::coefficients(quantreg::rq(y ~ 0 + x, tau = tau))
  }
  betaold = beta_all = beta = betaint


  if (is.null(f0)) {
    f0 = rep(1, length(y))
  }


  ul = x %*% beta

  ghat = ul
  iteration = 0
  w = 1
  diff2a = 1
  removed_ind = c()

  kk = 0 #maybe remove

  while (max(abs(diff2a)) > 0.001 & iteration < max_it) {
    iteration = iteration + 1

    new_removed_ind <- which(abs(beta) <= cutoff)
    if (length(new_removed_ind) != 0 &&
        length(c(removed_ind, new_removed_ind)) < -1 + length(beta_all)) {
      #update removed_ind
      if (is.null(removed_ind)) {
        removed_ind = new_removed_ind
      } else{
        removed_ind = c(removed_ind,
                        (1:length(beta_all))[-removed_ind][new_removed_ind])
      }
      #update X and beta
      if (!is.null(removed_ind)) {
        x = x_all[, -removed_ind]
        beta_all[removed_ind] = 0
        beta = beta_all[-removed_ind]
        beta_all[-removed_ind] = beta
      } else{
        x = x_all
        beta = beta_all
        beta_all = beta
      }
    }

    resid <- NULL
    R <- sumuv <- matrix(0, max(nk), max(nk))

    resid = 0.5 - (y - ghat < 0)
    ######## used to calculate the correlation matrix
    sd_resid = resid / sqrt(tau - tau ^ 2) ##########standard the resid

    if (worktype == "CS") {
      Ns = sum(nk * (nk - 1))
    } else if (worktype == "AR") {
      Ns = sum(nk - 1)
      Ns1 = sum(nk - 2)
    } else if (worktype == "Tri") {
      Ns = sum(nk - 1)
      Ns1 = sum(nk - 2)
    } else{
      Ns = nsub
    }

    sumrd <- sumrd1 <- 0
    R <- sumjk <- matrix(0, max(nk), max(nk))

    for (i in 1:nsub) {
      if (nk[i] >= 2) {
        za = sd_resid[(cn[i] + 1):cn[i + 1]]
        if (worktype == "CS") {
          for (j in 1:nk[i]) {
            for (k in 1:nk[i]) {
              if (k != j)
                sumrd = sumrd + za[j] * za[k]
            }
          }
        } else if (worktype == "AR") {
          for (j in 2:nk[i]) {
            sumrd = sumrd + za[j] * za[j - 1]
            if (j >= 3)
              sumrd1 = sumrd1 + za[j] * za[j - 2]
          }
        } else if (worktype == "Tri") {
          for (j in 2:nk[i])
            sumrd = sumrd + za[j] * za[j - 1]
        } else{
          for (j in 1:nk[i]) {
            for (k in 1:nk[i]) {
              R[j, k] = R[j, k] + za[j] * za[k]
              sumjk[j, k] = sumjk[j, k] + 1
            }
          }
        }
      }
    }

    if (worktype == "CS") {
      af = sumrd / Ns
      R = diag(max(nk))
      for (i in 1:max(nk)) {
        for (j in 1:max(nk)) {
          if (i != j)
            R[i, j] = af
        }
      }
    } else if (worktype == "AR") {
      #af=(sumrd/Ns+sqrt(max(sumrd1,0))/Ns1)/2;R=diag(max(nk))
      af = sumrd / Ns
      R = diag(max(nk))
      for (i in 1:max(nk)) {
        for (j in 1:max(nk)) {
          R[i, j] = af ^ (abs(i - j))
        }
      }
    } else if (worktype == "Tri") {
      af = (sumrd / Ns + sqrt(sumrd1) / Ns1) / 2

      R = diag(max(nk))
      for (i in 1:max(nk)) {
        for (j in 1:max(nk)) {
          if (abs(i - j) <= 1)
            R[i, j] = af ^ (abs(i - j))
        }
      }
    } else if (worktype == "Ind") {
      R = diag(max(nk))
    } else{
      R = R / (nsub - nx)

      temp <- sqrt(diag(R))
      R <- t(t(R / temp) / temp)

    }

    ############induce smoothing ##################
    U2 = 0
    H2 = 0

    r2 = sqrt(diag(x %*% t(x)))
    mu2 = x %*% beta
    hh = y - mu2
    fr = (pnorm(nsub ^ 0.5 * (y - mu2) / r2, 0, 1) - rep(1 - tau, N))


    for (i in 1:nsub) {
      a_ind = (cn[i] + 1):cn[i + 1]
      Ra = R[1:nk[i], 1:nk[i]]

      fra = fr[a_ind]
      r2a = r2[a_ind]
      ha = hh[a_ind]
      if (nk[i] == 1) {
        Ga = as.vector(f0[a_ind])
        xa = matrix(x[a_ind, ], nrow = 1)

        Ha2 = as.vector(dnorm(nsub ^ 0.5 * ha / r2a, 0, 1) / r2a)
      } else{
        Ga = diag(as.vector(f0[a_ind]))
        xa = x[a_ind, ]

        Ha2 = diag(as.vector(dnorm(nsub ^ 0.5 * ha / r2a, 0, 1) / r2a))
      }
      U2 = U2 + t(xa) %*% Ga %*% geninv(Ra) %*% fra
      H2 = H2 + nsub ^ 0.5 * t(xa) %*% Ga %*% geninv(Ra) %*% Ha2 %*% xa
    }



    if (lambda > 0) {
      #penalty
      eps <- 10 ^ -6
      #scad
      pe <-
        as.vector(pp_scad_sim(abs(as.vector(beta)), lambda) / (abs(as.vector(beta)) +
                                                                 eps))
      if (length(pe) != 1) {
        sigma <- nsub * diag(pe)
      } else{
        sigma <- nsub * pe
      }
      diff2 = geninv(H2 + sigma) %*% (U2 - sigma %*% beta)
    } else{
      diff2 = geninv(H2) %*% U2
    }
    beta = beta + w * diff2
    ghat = x %*% beta
    mcl = mean(check_loss(y - ghat, tau))

    if (max(abs(diff2)) > 100)
      break

    diff2a = max(abs(w * diff2))
    w = w / 2
  }

  converge = TRUE
  if (kk >= 3 || iteration >= max_it || max(abs(diff2)) > 100) {
    converge = FALSE
  }

  if (!is.null(removed_ind)) {
    beta_out = beta_all
    beta_out[removed_ind] = 0
    beta_out[-removed_ind] = beta
    X_selected = setdiff((1:length(beta_all)), paste(removed_ind))
  } else{
    beta_out = beta_all
    X_selected = (1:length(beta_all))
  }

  #calculate hbic
  hbic = log(mcl * N) + (log(nsub) / (2 * nsub)) * log(log(NCOL(x))) * length(X_selected)


  qpgee.obj = list(
    beta = beta_out,
    fitted = ghat,
    R = R,
    X_selected = X_selected,
    mcl = mcl,
    hbic = hbic,
    converge = converge
  )
  class(qpgee.obj) <- "qpgee"
  return(qpgee.obj)

}

#' Predict method for qpgee model objects
#'
#' This function makes predictions from a "qpgee" model object. When `newdata` is not provided,
#' it returns predictions using the original data the model was fitted on. If `newdata` is supplied
#' (through `...`), it uses this new data for prediction.
#'
#' @param object A "qpgee" model object.
#' @param ... Additional arguments to the function. Can include `newdata`, a dataframe
#'    containing the new data to predict on. The structure of `newdata` should match that of the data
#'    the model was originally fitted with, specifically in terms of the variables it contains.
#'    Additional arguments are ignored.
#'
#' @return If `newdata` is not supplied, returns a vector of predictions based on the fitted values
#'    and handling of NAs specified in the model object. If `newdata` is supplied, returns a vector
#'    of predictions for the new data.
#' @examples
#' # Example usage:
#' sim_data <- generateData(n_sub = 100, n_obs = rep(10, 100),  p = 100,
#'                          beta0 = rep(1,7), rho = 0.6, type = "ar",
#'                           dis = "normal", ka = 1)
#'
#' X=sim_data$X
#' y=sim_data$y
#'
#' #fit qpgee
#' qpgee.fit = qpgee(X,y,tau=0.5,nk=rep(10, 100))
#' predict(qpgee.fit)
#'
#' @export
predict.qpgee <- function (object, ...)
{
  args <- list(...)
  if (!"newdata" %in% names(args)) {
    return(stats::napredict(object$na.action, object$fitted))
  } else {
    X <- args[["newdata"]]
  }
  pred <- drop(X %*% object$beta)
  pred
}
