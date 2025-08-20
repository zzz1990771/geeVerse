
library(geeVerse)
#load libraries
rm(list=ls())
library(splines)
library(SparseM)
library(quantreg)
library(mvtnorm)
library(MASS)
#library(np)
library(foreach)
library(doParallel)
library(glmnet)
source("./R/plrq.R")
source("./R/qpgeeSim_new.R")
source("./R/qsicm.R")
source("./R/utils.R")



## Monte Carlo Simulation Script for Table 2 -----------------------------
n_sim = 1 # number of simulation replications
res_list = NULL # initialize a list to store all the fitted models
res_list_new = NULL
# start simulation
for (m in 1:n_sim){
  #generate data
  nsub = 100
  set.seed(2024)
  sim_data = generate_data(nsub = nsub, nobs = rep(10, nsub), p = 200,
                           beta0 = c(rep(1, 7), rep(0, 193)),
                           rho = 0.6, corstr = "exchangeable")
  # fit qpgee on different correlation structure for different tau
  correlation_vec = c("independence","exchangeable","AR1")
  tau_vec = c(0.1,0.5,0.9)
  for (correlation in correlation_vec){
    for(tau in tau_vec){
      QPGEE_results = qpgee(y ~ 0 + . -id, data = sim_data, id = sim_data$id,
                            corstr = correlation, tau = tau, ncore = 10,
                            control = qpgeeControl(shrinkCutoff = 10^-3))

      res_list[[correlation]][[as.character(tau)]][[m]] = QPGEE_results

      y = sim_data$y
      x = sim_data[,-c(1,2)]
      conversion_map <- c("independence" = "Ind", "exchangeable" = "CS", "AR1" = "AR")
      correlation_new <- conversion_map[correlation]
      aa2_scad <- qpgeePlsim_hbic_new(X=NULL,y=y,Z=x,
                                      tau = tau,
                                      method = "HBIC",
                                      ncore = 10,
                                      nobs = rep(10, nsub),
                                      corstr=correlation_new,
                                      #intercept = FALSE,
                                      f0 = rep(1,N),#NULL,
                                      max_it = 500,
                                      betaint=NULL,
                                      lambda=0.0,     # penalty for linear coefficients ##for beta (single-index)
                                      lambda_thta=0.0, # penalty for spline coefficients
                                      lambda_gamma=NULL ,# penalty for linear coefficients
                                      fast = TRUE,
                                      check_converge = FALSE)#cutoff = .001)
      aa2_scad$beta <- aa2_scad$gamma
      res_list_new[[correlation]][[as.character(tau)]][[m]] = aa2_scad

    }
  }
}
saveRDS(res_list,"res_list_old.RDS")
saveRDS(res_list_new,"res_list_new.RDS")
#Report results for independence working correlation structure with tau = 0.5 from Table 2;
#Note: users can update parameters for other working correlation structures and tau results
compile_result(res_list[["independence"]][[as.character(0.9)]],
               beta0= c(rep(1,7),rep(0,193)))
compile_result(res_list_new[["independence"]][[as.character(0.9)]],
               beta0= c(rep(1,7),rep(0,193)))

