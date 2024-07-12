#supporting functions

#correlated error function
Siga_cov<-function(rho,correlation,nt){
  sigma=matrix(0,nt,nt)
  if (correlation=="cs"){
    sigma=(1-rho)*diag(nt)+rho*matrix(1,nt,nt)
  }else if (correlation=="AR1"){
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


#compile results
compile_result <- function(PQGEE_results,beta0,cutoff = 0.1){
  p = length(beta0)
  n_sim = length(PQGEE_results)
  results_table <- matrix(NA,ncol=5,nrow=2)
  colnames(results_table) <- c("percentage","TP","FP","MSE","MAD")
  result_structure = matrix(NA,nrow=n_sim,ncol=5)
  for(sim in 1:n_sim){
    selected <- rep(FALSE,p)
    selected <- abs(PQGEE_results[[structure]][[sim]]$beta)>cutoff
    selected_true <- abs(beta0)>cutoff
    correct <- FALSE
    TP <- sum(selected & selected_true)
    FP <- sum(selected & !selected_true)
    if((TP-FP)==sum(abs(beta0)>0)){correct = TRUE}
    beta_est <- rep(0,p)
    beta_est <- PQGEE_results[[structure]][[sim]]$beta
    MSE <- sum((beta_est-beta)^2)
    MAD <- sum(abs(beta_est-beta))
    result_structure[sim,] <- c(correct,TP,FP,MSE,MAD)
  }
  results_table <- rbind(apply(result_structure,2,mean),
                         apply(result_structure,2,sd))

  return(cbind(results_table[1:2,],results_table[3:4,],results_table[5:6,]) )
}



#compile results
compile_result_hd <- function(PQGEE_results,SIS_results,beta0,cutoff = 0.1){
  p = length(beta0)
  n_sim = length(PQGEE_results)
  results_table <- matrix(NA,ncol=5,nrow=2)
  colnames(results_table) <- c("percentage","TP","FP","MSE","MAD")
  result_structure = matrix(NA,nrow=n_sim,ncol=5)
  for(sim in 1:n_sim){
    #SIS results
    SIS_X <- SIS_results[[sim]]
    selected <- rep(FALSE,p)
    selected[SIS_X] <- abs(PQGEE_results[[structure]][[sim]]$beta)>cutoff
    selected_true <- abs(beta0)>cutoff
    correct <- FALSE
    TP <- sum(selected & selected_true)
    FP <- sum(selected & !selected_true)
    if((TP-FP)==sum(abs(beta0)>0)){correct = TRUE}
    beta_est <- rep(0,p)
    beta_est[SIS_X] <- PQGEE_results[[structure]][[sim]]$beta
    MSE <- sum((beta_est-beta)^2)
    MAD <- sum(abs(beta_est-beta))
    result_structure[sim,] <- c(correct,TP,FP,MSE,MAD)
  }
  results_table <- rbind(apply(result_structure,2,mean),
                         apply(result_structure,2,sd))

  return(cbind(results_table[1:2,],results_table[3:4,],results_table[5:6,]) )

}


#scad penalty function
pp_scad <- function(z,lambda_si,a=3.7) {
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
