#supporting functions

#correlated error function
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


#compile results
compile_result <- function(PQGEE_results,n_sim,cutoff = 0.1,beta,beta0){
  results_table <- matrix(NA,ncol=5,nrow=6)
  colnames(results_table) <- c("percentage","TP","FP","MSE","MAD")
  i=1
  for(structure in c("Ind","CS","AR")){#structure="Ind" #sim=1
    result_structure = matrix(NA,nrow=n_sim,ncol=5)
    for(sim in 1:n_sim){
      selected <- abs(PQGEE_results[[structure]][[sim]]$beta)>cutoff
      selected_true <- abs(beta)>cutoff
      correct <- FALSE
      TP <- sum(selected & selected_true)
      FP <- sum(selected & !selected_true)
      if((TP-FP)==length(beta0)){correct = TRUE}
      MSE <- sum((PQGEE_results[[structure]][[sim]]$beta-beta)^2)
      MAD <- sum(abs(PQGEE_results[[structure]][[sim]]$beta-beta))
      result_structure[sim,] <- c(correct,TP,FP,MSE,MAD)
    }
    results_table[c(i,i+1),] <- rbind(apply(result_structure,2,mean),
                                      apply(result_structure,2,sd))
    i=i+2
  }

  return(cbind(results_table[1:2,],results_table[3:4,],results_table[5:6,]))

}








compile_result_one <- function(PQGEE_results,cutoff = 0.1,beta,beta0,structures=c("Ind","CS","AR")){
  n_sim=1
  results_table <- matrix(NA,ncol=7,nrow=6)
  colnames(results_table) <- c("F1 score","TP","FP","TN", "FN","MSE","MAD")
  i=1
  for(structure in structures ){#structure="Ind" #sim=1
    result_structure = matrix(NA,nrow=n_sim,ncol=7)
    for(sim in 1:n_sim){
      selected <- abs(PQGEE_results[[structure]][[sim]]$beta)>cutoff
      selected_true <- abs(beta)>cutoff
      correct <- FALSE
      TP <- sum(selected & selected_true)
      FP <- sum(selected & !selected_true)
      if((TP-FP)==length(beta0)){correct = TRUE}
      
      FN <-sum(!selected & selected_true)
      TN <-sum(!selected & !selected_true)
      #if((TP-FP)==length(beta0)){correct = TRUE}
      Precision = TP/(TP+FP)
      Recall = TP/(TP+FN)
      F1 = 2*(Precision*Recall)/(Precision+Recall)
      MSE <- sum((PQGEE_results[[structure]][[sim]]$beta-beta)^2)
      MAD <- sum(abs(PQGEE_results[[structure]][[sim]]$beta-beta))
      result_structure[sim,] <-  c(F1, TP,FP,TN,FN,MSE,MAD)
      #c(correct,TP,FP,MSE,MAD)
    }
    results_table[c(i),] <- result_structure
    i=i+2
  }
  
  return(cbind(results_table[1:2,],results_table[3:4,],results_table[5:6,]))
  
}
 


#compile results
compile_result_hd <- function(PQGEE_results,SIS_results,n_sim,cutoff = 0.1,beta,beta0,p){
  results_table <- matrix(NA,ncol=5,nrow=6)
  colnames(results_table) <- c("percentage","TP","FP","MSE","MAD")
  i=1
  for(structure in c("Ind","CS","AR")){
    result_structure = matrix(NA,nrow=n_sim,ncol=5)
    for(sim in 1:n_sim){
      #SIS results
      SIS_X <- SIS_results[[sim]]
      selected <- rep(FALSE,p)
      selected[SIS_X] <- abs(PQGEE_results[[structure]][[sim]]$beta)>cutoff
      selected_true <- abs(beta)>cutoff
      correct <- FALSE
      TP <- sum(selected & selected_true)
      FP <- sum(selected & !selected_true)
      if((TP-FP)==length(beta0)){correct = TRUE}
      beta_est <- rep(0,p)
      beta_est[SIS_X] <- PQGEE_results[[structure]][[sim]]$beta
      MSE <- sum((beta_est-beta)^2)
      MAD <- sum(abs(beta_est-beta))
      result_structure[sim,] <- c(correct,TP,FP,MSE,MAD)
    }
    results_table[c(i,i+1),] <- rbind(apply(result_structure,2,mean),
                                      apply(result_structure,2,sd))
    i=i+2
  }
  return(cbind(results_table[1:2,],results_table[3:4,],results_table[5:6,]) )

}

compile_result_one_sim <- function(PQGEE_results,SIS_results,n_sim,cutoff = 0.1,beta,beta0,p){
  results_table <- matrix(NA,ncol=7,nrow=6)
  colnames(results_table) <- c("F1 score","TP","FP","TN", "FN","MSE","MAD")
  i=1
  for(structure in c("Ind","CS","AR")){
    result_structure = matrix(NA,nrow=n_sim,ncol=7)
    for(sim in 1:n_sim){
      # #SIS results
      # #SIS_X <- SIS_results[[sim]]
      # selected <- rep(FALSE,p)
      # selected[SIS_X] <- abs(PQGEE_results[[structure]][[sim]]$beta)>cutoff
      # selected_true <- abs(beta)>cutoff
      
      selected <- abs(PQGEE_results[[structure]][[sim]]$beta)>cutoff
      selected_true <- abs(beta)>cutoff
      #correct <- FALSE
      # correct <- FALSE
      TP <- sum(selected & selected_true)
      FP <- sum(selected & !selected_true)
      FN <-sum(!selected & selected_true)
      TN <-sum(!selected & !selected_true)
      #if((TP-FP)==length(beta0)){correct = TRUE}
      Precision = TP/(TP+FP)
      Recall = TP/(TP+FN)
      F1 = 2*(Precision*Recall)/(Precision+Recall)
      MSE <- sum((PQGEE_results[[structure]][[sim]]$beta-beta)^2)
      MAD <- sum(abs(PQGEE_results[[structure]][[sim]]$beta-beta))
      result_structure[sim,] <- c(F1, TP,FP,TN,FN,MSE,MAD)
    }
    results_table[i,] <- apply(result_structure,2,mean)
    i=i+2
  }
  output=cbind(results_table[1,],results_table[3,],results_table[5,])
  colnames(output)<-c("Ind","CS","AR")
  output= round(output,3)
  return(output) 
  
}



compile_result_demonstration_CCI<- function(PQGEE_results){
 #PQGEE_results=PQGEE_resultsMedian
   #results_table <- matrix(NA,ncol=7,nrow=6)
  
  result_structure = matrix(NA,nrow=3,ncol=2)
  colnames(result_structure) <- c("MCL","# selected")
  WCM =c("Ind","CS","AR")
  for(structure in WCM){
    
    struct_num=match(structure, WCM)
     result_structure[struct_num,] <- c(PQGEE_results[[struct_num]]$mcl,
                                        length(PQGEE_results[[struct_num]]$X_selected))
    # results_table <- apply(result_structure,2,mean)
}

  #output=cbind(results_table[1,],results_table[3,],results_table[5,])
  #colnames(output)<-structure
  output= round(result_structure,4)
  return(output) 
  
}




compile_result_hd_one_sim <- function(PQGEE_results,SIS_results,n_sim,cutoff = 0.1,beta,beta0,p){
  results_table <- matrix(NA,ncol=7,nrow=6)
  colnames(results_table) <- c("F1 score","TP","FP","TN", "FN","MSE","MAD")
  i=1
  for(structure in c("Ind","CS","AR")){
    result_structure = matrix(NA,nrow=n_sim,ncol=7)
    for(sim in 1:n_sim){
      #SIS results
      SIS_X <- SIS_results[[sim]]
      selected <- rep(FALSE,p)
      selected[SIS_X] <- abs(PQGEE_results[[structure]][[sim]]$beta)>cutoff
      selected_true <- abs(beta)>cutoff
     # correct <- FALSE
      TP <- sum(selected & selected_true)
      FP <- sum(selected & !selected_true)
      FN <-sum(!selected & selected_true)
      TN <-sum(!selected & !selected_true)
      #if((TP-FP)==length(beta0)){correct = TRUE}
      Precision = TP/(TP+FP)
      Recall = TP/(TP+FN)
      F1 = 2*(Precision*Recall)/(Precision+Recall)
      beta_est <- rep(0,p)
      beta_est[SIS_X] <- PQGEE_results[[structure]][[sim]]$beta
      MSE <- sum((beta_est-beta)^2)
      MAD <- sum(abs(beta_est-beta))
      result_structure[sim,] <- c(F1, TP,FP,TN,FN,MSE,MAD)
    }
    results_table[i,] <- apply(result_structure,2,mean)
    i=i+2
  }
  output=cbind(results_table[1,],results_table[3,],results_table[5,])
  colnames(output)<-c("Ind","CS","AR")
  output= round(output,3)
  return(output) 
  
}


compile_result_rea_data <- function(PQGEE_results,SIS_results,n_sim,cutoff = 0.1,beta,beta0,p){
  results_table <- matrix(NA,ncol=7,nrow=6)
  colnames(results_table) <- c("F1 score","TP","FP","TN", "FN","MSE","MAD")
  i=1
  for(structure in c("Ind","CS","AR")){
    result_structure = matrix(NA,nrow=n_sim,ncol=7)
    for(sim in 1:n_sim){
      #SIS results
      SIS_X <- SIS_results[[sim]]
      selected <- rep(FALSE,p)
      selected[SIS_X] <- abs(PQGEE_results[[structure]][[sim]]$beta)>cutoff

      beta_est[SIS_X] <- PQGEE_results[[structure]][[sim]]$beta
     # MSE <- sum((beta_est-beta)^2)
     # MAD <- sum(abs(beta_est-beta))
      result_structure[sim,] <- c(F1, TP,FP,TN,FN,MSE,MAD)
    }
    results_table[i,] <- apply(result_structure,2,mean)
    i=i+2
  }
  output=cbind(results_table[1,],results_table[3,],results_table[5,])
  colnames(output)<-c("Ind","CS","AR")
  output= round(output,3)
  return(output) 
  
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
