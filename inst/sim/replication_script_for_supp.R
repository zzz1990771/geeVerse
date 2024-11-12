# This replication script includes all supplementary simulation codes
# for manuscript `geeVerse: Ultra-high Dimensional Heterogeneous Data Analysis
# with Generalized Estimating Equations`. The results in the Monte Carlo
# simulation study may be slightly different due to randomness.


# Before you run the script, make sure you have the geeVerse package installed.
# It can be installed from a source package provided or directly from CRAN.
library(geeVerse)
library(tictoc)

tic()

## Example 1 Monte Carlo Simulation Script-----------------------------------------
n_sim = 100 # number of simulation replications
res_list = NULL # initialize a list to store all the fitted models
# start simulation
for (m in 1:n_sim){
  #generate data
  nsub = 100
  sim_data = generateData(nsub= nsub, nobs = rep(10, nsub),  p = 200,
        beta0 = c(rep(1,7),rep(0,193)), rho = 0.6, correlation = "exchangeable")
  y = sim_data$y
  x = sim_data$X
  # fit qpgee on different correlation structure for different tau
  correlation_vec = c("independence","exchangeable","AR1")
  tau_vec = c(0.1,0.5,0.9)
  for (correlation in correlation_vec){
    for(tau in tau_vec){
      QPGEE_results = qpgee(x, y, nobs=rep(10, nsub), correlation = correlation,
                            tau = tau, max_it=100, method = "HBIC",ncore = 10)
      res_list[[correlation]][[as.character(tau)]][[m]] = QPGEE_results
    }
  }
}
#Report results for independence working correlation structure with tau = 0.5 from Table 2;
#Note: users can update parameters for other working correlation structures and tau results
compile_result(res_list[["independence"]][[as.character(0.5)]],
               beta0= c(rep(1,7),rep(0,193)))

## Example 2 -----------------------------------------


## Monte Carlo Simulation  Script -----------------------------
n_sim = 100 # number of simulation replications
res_list = NULL # initialize a list to store all the fitted models
# start simulation
for (m in 1:n_sim){
  #generate data
  nsub = 100
  sim_data = generateData(nsub= nsub, nobs = rep(10, nsub),  p = 200,
                          beta0 = c(rep(1,7),rep(0,193)), rho = 0.6, correlation = "indepedence")
  y = sim_data$y
  x = sim_data$X
  # fit qpgee on different correlation structure for different tau
  correlation_vec = c("independence")
  tau_vec = c(0.1,0.5,0.9)
  for (correlation in correlation_vec){
    for(tau in tau_vec){
      QPGEE_results = qpgee(x, y, nobs=rep(1, 1000), correlation = correlation,
                            tau = tau, max_it=100, method = "HBIC",ncore = 10)
      res_list[[correlation]][[as.character(tau)]][[m]] = QPGEE_results
    }
  }
}
#Report results for independence working correlation structure with tau = 0.5 from Table 3;
#Note: users can update parameters for other working correlation structures and tau results
compile_result(res_list[["independence"]][[as.character(0.5)]],
               beta0= c(rep(1,7),rep(0,193)))

## Example 3 -----------------------------------------

## Monte Carlo Simulation Script ---------------------------------
n_sim = 100 # number of simulation replications
res_list = NULL # initialize a list to store all the fitted models
# start simulation
for (m in 1:n_sim){
  #generate data
  nsub = 100
  nobs = c(rep(10,nsub/2),rep(5,nsub/2))
  sim_data = generateData(nsub= 100, nobs = nobs, p = 200,
                          beta0 = c(rep(1,7),rep(0,193)), rho = 0.6, correlation = "exchangeable")
  y = sim_data$y
  x = sim_data$X
  # fit qpgee on different correlation structure for different tau
  correlation_vec = c("independence","exchangeable","AR1")
  tau_vec = c(0.1,0.5,0.9)
  for (correlation in correlation_vec){
    for(tau in tau_vec){
      QPGEE_results = qpgee(x, y, nobs=nobs, correlation = correlation,
                            tau = tau, max_it=100, method = "HBIC",ncore = 10)
      res_list[[correlation]][[as.character(tau)]][[m]] = QPGEE_results
    }
  }
}
#Report results for independence working correlation structure with tau = 0.5 from Table 4;
#Note: users can update parameters for other working correlation structures and tau results
compile_result(res_list[["independence"]][[as.character(0.5)]],
               beta0= c(rep(1,7),rep(0,193)))


# Extended Monte Carlo Simulation -------------------------------------------------------------
# This subsection includes codes for table 5 and table 6.
# Note we use 2 replications for computational efficiency instead of 100 replications
# as in the paper, so our results will likely be different than the 100 replication simulation based tables in the paper.



## Monte Carlo Simulation Script for Table 5 ----------------------------------
n_sim = 100 # number of simulation replications
res_list = NULL # initialize a list to store all the fitted models
# start simulation
for (m in 1:n_sim){
  #generate data
  nsub = 100
  sim_data = generateData(nsub= nsub, nobs = rep(10, nsub),  p = 200,
                          beta0 = c(rep(1,7),rep(0,193)), rho = 0.6, correlation = "exchangeable")
  y = sim_data$y
  x = sim_data$X
  # fit qpgee on different correlation structure for different tau
  correlation_vec = c("independence","exchangeable","AR1")
  tau_vec = c(0.1,0.5,0.9)
  for (correlation in correlation_vec){
    for(tau in tau_vec){
      QPGEE_results = qpgee(x, y, nobs=rep(10, nsub), correlation = correlation,
                            tau = tau, max_it=100, method = "HBIC",ncore = 10)
      res_list[[correlation]][[as.character(tau)]][[m]] = QPGEE_results
    }
  }
}
#Report results for independence working correlation structure with tau = 0.5 from Table 5;
#Note: users can update parameters for other working correlation structures
compile_result(res_list[["independence"]][[as.character(0.5)]],
               beta0= c(rep(1,7),rep(0,193)))


## Monte Carlo Simulation Script for Table 6 ---------------------------
n_sim = 1 # number of simulation replications
res_list = NULL # initialize a list to store all the fitted models
# start simulation
for (m in 1:n_sim){
  #generate data
  nsub = 100
  set.seed(2024)
  sim_data = generateData(nsub = nsub, nobs = rep(10, nsub), p = 1000,
                          beta0 = c(rep(1,7),rep(0,993)), rho = 0.6, correlation = "exchangeable")
  y = sim_data$y
  x = sim_data$X
  # fit qpgee on different correlation structure for different tau
  correlation_vec = c("independence","exchangeable","AR1")
  tau_vec = c(0.1,0.5,0.9)
  for (correlation in correlation_vec){
    for(tau in tau_vec){
      x_sis_ind = SIS::SIS(x,y,nsis=200)$sis.ix0
      xsis = x[,x_sis_ind]
      QPGEE_results = qpgee(xsis, y, nobs=rep(10, nsub), correlation = correlation,
                            tau = tau, max_it=100, method = "HBIC",ncore = 10)
      full_beta = rep(0,1000)
      full_beta[x_sis_ind] = QPGEE_results$beta
      QPGEE_results$beta = full_beta
      res_list[[correlation]][[as.character(tau)]][[m]] = QPGEE_results
    }
  }
}
#Report results for independence working correlation structure with tau = 0.5 from Table 6;
#Note: users can update parameters for other working correlation structures
compile_result(res_list[["exchangeable"]][[as.character(0.9)]],
               beta0= c(rep(1,7),rep(0,993)))

## Computational Time Comparison Script, Table 7 --------------------------

# Note that we only showcase p = 200 here.
# For p = 2000, just change p from 200 to 2000.
p = 200
# Initialize some vectors to store computation time
# No palatalization is used for all three methods tested
qpgee_tvec = c()
pgee_own_tvec = c()
pgee_pkg_tvec = c()

# start simulation
n_sim = 1
for (i in 1:n_sim) {
  #generate data
  nsub = 100
  nobs = rep(10, nsub)
  id = rep(1:nsub,each = nobs)
  sim_data = generateData(nsub= nsub, nobs = rep(10, nsub),  p = p,
                          beta0 = c(rep(1,7),rep(0,193)), rho = 0.6, correlation = "exchangeable")
  y = sim_data$y
  x = sim_data$X

  #fit qpgee
  qpgee_time = system.time(qpgee_fit <-
                             qpgee(x,y,lambda = 0.1, tau=0.5,nobs=nobs,cutoff = 10^-3))[3]
  qpgee_tvec = c(qpgee_tvec,qpgee_time)

  #fit own pgee
  data = data.frame(x,y,id)
  pgee_own_time = system.time(PGEE_own_fit <- PGEE("y ~.-id-1",id = id,data = data,
                                                   corstr = "exchangeable",lambda=0.1))[3]
  pgee_own_tvec = c(pgee_own_tvec,pgee_own_time)

  #fit pgee from PGEE package
  pgee_pkg_time = system.time(PGEE_pkg <- PGEE::PGEE("y ~.-id-1",id = id,data = data,
                                                     corstr = "exchangeable",lambda=0.1))[3]
  pgee_pkg_tvec = c(pgee_pkg_tvec,pgee_pkg_time)

}
#a tool function to report performance
report_sum <- function(x){
  c(mean(x),min(x),max(x),median(x))
}
#report results for table 7
report_sum(qpgee_tvec)
report_sum(pgee_own_tvec)
report_sum(pgee_pkg_tvec)


## Monte Carlo Simulation Script for Table 8 ---------------------------
# This subsection includes codes for Table 8, p = 50
n_sim = 1 # number of simulation replications
res_list = NULL # initialize a list to store all the fitted models
# start simulation
for (m in 1:n_sim){
  #generate data
  nsub = 1000
  nobs = rep(5,1000)
  p = 50
  sim_data <- generateData(nsub = nsub, nobs = nobs, p = p,
                           beta0 = c(rep(1,9),rep(0,41)), rho = 0.6, correlation = "exchangeable",
                           ka = 0.5, SNPs = simuGene[,1:25])
  y = sim_data$y
  x = sim_data$X
  #for PGEE, reorganize the data
  id = rep(1:nsub, times = nobs)
  data = data.frame(x,y,id)
  # fit qpgee on different correlation structure for different tau
  correlation_vec = c("independence","exchangeable","AR1")
  #correlation_vec = c("exchangeable","AR1")
  for (correlation in correlation_vec){
    # qpgee tau = 0.5
    tau = 0.5
    QPGEE_results = qpgee(x, y, nobs=nobs, correlation = correlation,
                          tau = tau, max_it=100, method = "HBIC",ncore = 18)
    res_list[[correlation]][["qpgee"]][[m]] = QPGEE_results
  }
}
#Report results for independence working correlation structure with tau = 0.5 from Table 8;
#Note: users can update parameters for other working correlation structures
compile_result(res_list[["exchangeable"]][["qpgee"]],
               beta0= c(rep(1,9),rep(0,41)))


# This subsection includes codes for Table 8, p = 1000
n_sim = 1 # number of simulation replications
res_list = NULL # initialize a list to store all the fitted models
# start simulation
for (m in 1:n_sim){
  #generate data
  nsub = 1000
  nobs = rep(5,1000)
  p = 1000
  sim_data <- generateData(nsub = nsub, nobs = nobs, p = p,
                           beta0 = c(rep(1,9),rep(0,991)), rho = 0.6, correlation = "exchangeable",
                           ka = 0.5, SNPs = simuGene)
  y = sim_data$y
  x = sim_data$X

  #SIS
  x_sis_ind = SIS::SIS(x,y,nsis=200)$sis.ix0
  xsis = x[,x_sis_ind]

  #for PGEE, reorganize the data
  id = rep(1:nsub, times = nobs)
  data = data.frame(xsis,y,id)
  # fit qpgee on different correlation structure for different tau
  correlation_vec = c("independence","exchangeable","AR1")
  for (correlation in correlation_vec){
    # fit qpgee with tau = 0.5
    tau = 0.5

    QPGEE_results = qpgee(xsis, y, nobs=nobs, correlation = correlation,
                          tau = tau, max_it=100, method = "HBIC",ncore = 10)
    # map screened beta back to full length
    full_beta = rep(0,p)
    full_beta[x_sis_ind] = QPGEE_results$beta
    QPGEE_results$beta = full_beta
    res_list[[correlation]][["qpgee"]][[m]] = QPGEE_results
  }
}
#Report results for independence working correlation structure with tau = 0.5 from Table 8;
#Note: users can update parameters for other working correlation structures
compile_result(res_list[["exchangeable"]][["qpgee"]],
               beta0= c(rep(1,9),rep(0,991)))

