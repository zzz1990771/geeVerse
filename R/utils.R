#scad penalty function
pp_scad_sim <- function(z,lambda_si,lambda_pl=0,a=3.7) {
  # sim model so dsi always =length(z)
  dsi=length(z)
  x <- matrix(abs(z), ncol=1)
  value <- matrix(c(rep(0,dim(x)[[1]])),ncol=1)
  #penalize separately here
  for(i in 1:dim(x)[[1]]){
    lambda<-ifelse(i<=dsi,lambda_si,lambda_pl)
    # print(sprintf("coeff is %g and lambda is %g",i,lambda))

    if(x[i] < lambda){value[i] <- lambda}
    else if (x[i] < a*lambda){
      value[i] <- (a*lambda-x[i])/(a-1)}
    else{
      value[i] <- 0}
  }
  return(value)
}


#correlated error function
#' Generate Covariance Matrix
#'
#' This function generates a covariance matrix based on the specified correlation structure.
#' The function supports "compound symmetry" (cs) and "autoregressive" (ar) correlation structures,
#' as well as an identity matrix as the default option when neither "cs" nor "ar" is specified.
#'
#' @param rho Numeric, the correlation coefficient used for generating the covariance matrix.
#'        For "cs", it represents the common correlation between any two observations.
#'        For "ar", it represents the correlation between two consecutive observations,
#'        with the correlation decreasing for observations further apart.
#' @param type Character, specifies the type of correlation structure for the covariance matrix.
#'        Options are "cs" for compound symmetry, "ar" for autoregressive, and any other input
#'        will result in an identity matrix.
#' @param nt Integer, the dimension of the square covariance matrix (number of time points or observations).
#'
#' @return A square matrix of dimension `nt` representing the specified covariance structure.
#'
#' @export
Siga_cov<-function(rho,type,nt){
  sigma=matrix(0,nt,nt)
  if (type=="cs"){
    sigma=(1-rho)*diag(nt)+rho*matrix(1,nt,nt)
  }else if (type=="ar"){
    for (i in 1:nt)
      for (j in 1:nt)
        sigma[i,j]=rho^(abs(i-j))
  }else{
    sigma=diag(nt)
  }
  return(sigma)
}

#check loss
check_loss <- function(u,tau){
  u*(tau-(u<0))
}

#log space between 10^a to 10^b
logspace <- function( d1, d2, n=20) exp(log(10)*seq(d1, d2, length.out=n))


# #compile results
# compile_result <- function(PQGEE_results,n_sim,cutoff = 0.1){
#   results_table <- matrix(NA,ncol=5,nrow=6)
#   colnames(results_table) <- c("percentage","TP","FP","MSE","MAD")
#   i=1
#   for(structure in names(PQGEE_results)){
#     result_structure = matrix(NA,nrow=n_sim,ncol=5)
#     for(sim in 1:n_sim){
#       selected <- abs(PQGEE_results[[structure]][[sim]]$beta)>cutoff
#       selected_true <- abs(beta)>cutoff
#       correct <- FALSE
#       TP <- sum(selected & selected_true)
#       FP <- sum(selected & !selected_true)
#       if((TP-FP)==length(beta0)){correct = TRUE}
#       MSE <- sum((PQGEE_results[[structure]][[sim]]$beta-beta)^2)
#       MAD <- sum(abs(PQGEE_results[[structure]][[sim]]$beta-beta))
#       result_structure[sim,] <- c(correct,TP,FP,MSE,MAD)
#     }
#     results_table[c(i,i+1),] <- rbind(apply(result_structure,2,mean),
#                                       apply(result_structure,2,sd))
#     i=i+2
#   }
#   results_table <- format(round(results_table,digits=4),nsmall = 4)
#   results_table[,1:3] <- format(round(as.numeric(results_table[,1:3]),digits=2),nsmall = 2)
#   return(cbind(results_table[1:2,],results_table[3:4,],results_table[5:6,]))
#
# }

# #compile results
# compile_result_intercept <- function(PQGEE_results,n_sim,cutoff = 0.1){
#   results_table <- matrix(NA,ncol=5,nrow=6)
#   colnames(results_table) <- c("percentage","TP","FP","MSE","MAD")
#   i=1
#   for(structure in names(PQGEE_results)){
#     result_structure = matrix(NA,nrow=n_sim,ncol=5)
#     for(sim in 1:n_sim){
#       selected <- abs(PQGEE_results[[structure]][[sim]]$beta)[-1]>cutoff
#       selected_true <- abs(beta)>cutoff
#       correct <- FALSE
#       TP <- sum(selected & selected_true)
#       FP <- sum(selected & !selected_true)
#       if((TP-FP)==sum(beta!=0)){correct = TRUE}
#       MSE <- sum((PQGEE_results[[structure]][[sim]]$beta[-1]-beta)^2)
#       MAD <- sum(abs(PQGEE_results[[structure]][[sim]]$beta[-1]-beta))
#       result_structure[sim,] <- c(correct,TP,FP,MSE,MAD)
#     }
#     results_table[c(i,i+1),] <- rbind(apply(result_structure,2,mean),
#                                       apply(result_structure,2,sd))
#     i=i+2
#   }
#   results_table <- format(round(results_table,digits=4),nsmall = 4)
#   results_table[,1:3] <- format(round(as.numeric(results_table[,1:3]),digits=2),nsmall = 2)
#   return(cbind(results_table[1:2,],results_table[3:4,],results_table[5:6,]))
#
# }
#
# #compile results
# compile_result_hd <- function(PQGEE_results,SIS_results,n_sim,cutoff = 0.1){
#   results_table <- matrix(NA,ncol=5,nrow=6)
#   colnames(results_table) <- c("percentage","TP","FP","MSE","MAD")
#   i=1
#   for(structure in names(PQGEE_results)){
#     result_structure = matrix(NA,nrow=n_sim,ncol=5)
#     for(sim in 1:n_sim){
#       #SIS results
#       SIS_X <- SIS_results[[sim]]
#       selected <- rep(FALSE,p)
#       selected[SIS_X] <- abs(PQGEE_results[[structure]][[sim]]$beta)>cutoff
#       selected_true <- abs(beta)>cutoff
#       correct <- FALSE
#       TP <- sum(selected & selected_true)
#       FP <- sum(selected & !selected_true)
#       if((TP-FP)==length(beta0)){correct = TRUE}
#       beta_est <- rep(0,p)
#       beta_est[SIS_X] <- PQGEE_results[[structure]][[sim]]$beta
#       MSE <- sum((beta_est-beta)^2)
#       MAD <- sum(abs(beta_est-beta))
#       result_structure[sim,] <- c(correct,TP,FP,MSE,MAD)
#     }
#     results_table[c(i,i+1),] <- rbind(apply(result_structure,2,mean),
#                                       apply(result_structure,2,sd))
#     i=i+2
#   }
#   results_table <- format(round(results_table,digits=4),nsmall = 4)
#   results_table[,1:3] <- format(round(as.numeric(results_table[,1:3]),digits=2),nsmall = 2)
#   return(cbind(results_table[1:2,],results_table[3:4,],results_table[5:6,]) )
#
# }



cosf <- function(t, p) {
  (cos(t)) ^ (p - 2)
}

zeroDelt <- function(t, eps, pp) {
  stopifnot(t <= 1)
  given <- eps * stats::integrate(cosf, 0, pi / 2, p = pp)$value
  stats::integrate(cosf, 0, asin(t), p = pp)$value - given
}

computeDelta <- function(eps, p) {
  #compute delta value for Singular Value upper bound. Ref: Hochstenbach (2013)
  1 / stats::uniroot(zeroDelt, interval = c(0, 0.5), eps = eps, pp = p)$root
}
