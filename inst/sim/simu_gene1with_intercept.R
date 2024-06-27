#table 7
library(geeVerse)

# read simulated SNPs
data(simuDataGene)

# In this simulation, we have 1000 participants with 1000 variables where 9 are important variables.
dim(unique(t(X_snp)))
#settings
n=n_sub=1000
p=1000
beta0=rep(1,9)
p0=length(beta0)
beta = c(beta0,rep(0,p-p0))
n_obs<-rep(5,n_sub);

ka = 0.5
rho=0.5
type="ex"
dis="normal"
n_sim = 20
#Please note that n_sim=2 is just for quick demo


# Simulate the Data

# Repeat each row 5 times
X_snp_lon <- X_snp[rep(1:nrow(X_snp), each=5), ]


tau_list <- c(0.5)
final_results_sim <- matrix(NA,ncol=15,nrow=6)
start_time = Sys.time()
for(tau in tau_list){
  message(tau)
  PQGEE_results <- NULL
  PGEE_results <- NULL
  SIS_results <- NULL
  sim = 1
  while(sim <= n_sim){
    message(sim)
    try({
      #print(paste0("simulation number is ",sim))
      #generate errors for each subject
      e = NULL
      id<-NULL
      for (i in 1:n_sub){
        id<-c(id,rep(i,n_obs[i]))
        sigmai=Siga_cov(rho,type,n_obs[i])
        if (dis=="normal") ei=mvtnorm::rmvnorm(1, mean=rep(0, n_obs[i]), sigma=sigmai)
        if (dis=="t") ei=mvtnorm::rmvt(1, sigmai, df = 4, delta = rep(0, n_obs[i]))
        e=c(e,ei);
      }

      #generate y and X
      N=sum(n_obs)
      nk=n_obs
      cn = c(0, cumsum(n_obs))
      X_phone1 = matrix(rnorm(n*25),n*5,5)
      X_phone2 = matrix(rnorm(n*5*(p/2-5)),n*5,p/2-5)
      x= X = cbind(X_phone1,X_snp_lon,X_phone2)
      #x= X = cbind(X_phone1,X_phone2,X_snp_lon)
      y=X%*%beta+(1+ka*abs(X[,1]))*e


      ####implement SIS for dimension reduction
      #screening by SIS
      SIS_result <- SIS::SIS(x,y,iter=TRUE,nsis=50)
      selected_X <- SIS_result$ix0

      true_selected_vars<-which(abs(beta)>0)

      ##check to see how well SIS performs
      not_contained<-true_selected_vars[!(true_selected_vars %in% selected_X)]

      length(true_selected_vars)
      length(not_contained)

      SIS_results[[sim]] <- selected_X

      #keep screened X
      x=X[,selected_X]

      #non-longitudinal quantile regression as initial
      betaint=coefficients(quantreg::rq(y~x,tau = tau))

      #Apply proposed method with hbic tuning
      for(structure in c("Ind","CS","AR")){
        PQGEE_results[[structure]][[sim]] <-  qpgee_tune(x,y,tau=tau,method="HBIC",
                                                         betaint=betaint,nk=nk,cutoff = 0.05,
                                                         intercept = TRUE,
                                                         worktype=structure,ncore = 18)
        #remove intercept from beta for compiling the results
        PQGEE_results[[structure]][[sim]]$beta = PQGEE_results[[structure]][[sim]]$beta[-1]
      }

      #PGEE
      id = rep(1:n_sub, each = 5)
      data = data.frame(x,y,id)

      lambda.vec <- seq(0.1,1,0.1)
      #lambda.vec <- exp(seq(log(10),log(0.01),length.out = 30))

      ## single core version:
      #cv <- CVfit(formula = "y ~.-id", id = id, data = data,
      #            fold = 5, lambda.vec = lambda.vec, pindex = NULL, eps = 10^-6, maxiter = 30,
      #            tol = 10^-6)
      ## parallel version:
      cv <- CVfit(formula = "y ~.-id", id = id, data = data,
                  fold = 5, lambda.vec = lambda.vec, pindex = NULL, eps = 10^-6, maxiter = 30,
                  tol = 10^-6, ncore= 10)

      #Apply proposed method with hbic tuning
      for(structure in c("independence","exchangeable","AR-1")){
        res_obj <- NULL
        PGEE_fit = PGEE("y ~.-id",id = id, data = data,corstr = structure,lambda=cv$lam.opt)
        res_obj$beta <- PGEE_fit$coefficients[-1]
        PGEE_results[[structure]][[sim]] <-  res_obj
      }
      print("done!")
      sim = sim + 1
    })

  }

  #compile results
  compile_result_hd <- function(PQGEE_results,SIS_results,n_sim,cutoff = 0.1,beta,beta0,p){
    results_table <- matrix(NA,ncol=5,nrow=6)
    colnames(results_table) <- c("percentage","TP","FP","MSE","MAD")
    i=1
    for(structure in names(PQGEE_results)){
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

  # compile results ------------------------------------
  final_results_sim[1:2,] <- compile_result_hd(PQGEE_results,SIS_results,n_sim,cutoff = 0.1,beta=beta,beta0=beta0,p=p)
  final_results_sim[3:4,] <- compile_result_hd(PGEE_results,SIS_results,n_sim,cutoff = 0.1,beta=beta,beta0=beta0,p=p)
}


end_time = Sys.time()
total_time=end_time - start_time



# Print results
xtable::xtable(final_results_sim,
               digits=c(1,rep(c(2,2,2,4,4),3)))

knitr::kable(final_results_sim)


