#library
if(!require(mvtnorm)){install.packages("mvtnorm");library(mvtnorm)}else{library(mvtnorm)}
if(!require(quantreg)){install.packages("quantreg");library(quantreg)}else{library(quantreg)}
if(!require(MASS)){install.packages("MASS");library(MASS)}else{library(MASS)}
if(!require(foreach)){install.packages("foreach");library(foreach)}else{library(foreach)}
library(geeVerse)
#self-made functions
source("./inst/others/utils.R")



#settings
n=n_sub=400
p=7
beta0=rep(1,7)
p0=length(beta0)
beta = c(beta0,rep(0,p-p0))
n_obs<-rep(10,n_sub);

ka = 1
rho=0.6
type="ar"
dis="normal"
n_sim = 200
#Please note that n_sim=2 is just for quick demo


# Simulate the Data


tau_list <- c(0.1,0.5,0.9)
final_results_sim <- matrix(NA,ncol=15,nrow=6)
for(tau in tau_list){
  PQGEE_results <- NULL
  sim = 1
  while(sim <= n_sim){
    #generate errors for each subject
    e = NULL
    id<-NULL
    for (i in 1:n_sub){
      id<-c(id,rep(i,n_obs[i]))
      sigmai=Siga_cov(rho,type,n_obs[i])
      if (dis=="normal") ei=rmvnorm(1, mean=rep(0, n_obs[i]), sigma=sigmai)
      if (dis=="t") ei=rmvt(1, sigmai, df = 4, delta = rep(0, n_obs[i]))
      e=c(e,ei);
    }

    #generate y and X
    N=sum(n_obs)
    nk=n_obs
    cn = c(0, cumsum(n_obs))
    x=X=matrix(rnorm(N*p),N,p)
    y=X%*%beta+(1+ka*abs(X[,1]))*e
    try({
      #non-longitudinal quantile regression as initial
      betaint=coefficients(rq(y~0+X,tau = tau))

      #Apply proposed method with hbic tuning
      for(structure in c("Ind","CS","AR")){
        PQGEE_results[[structure]][[sim]] <-  qpgee(x,y,tau = tau, betaint = betaint,nk = nk,
                                                    worktype = structure,lambda = 0, cutoff = 0.01)
      }


      sim = sim +1
    })

  }

  # compile results ------------------------------------
  final_results_sim[c(2*which(tau==tau_list)-1,2*which(tau==tau_list)),] <- compile_result(PQGEE_results,n_sim,cutoff = 0.01,beta,beta0)

}




# Print results

xtable::xtable(final_results_sim,
               digits=c(1,rep(c(2,2,2,4,4),3)))

knitr::kable(final_results_sim)


#test geepack
# Fit a GEE model
#gee_model <- geeglm(y ~ x, id = id, family = gaussian, corstr = "exchangeable")

