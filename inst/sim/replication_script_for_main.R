# This replication script includes all demonstration codes
# for manuscript `geeVerse: Ultra-high Dimensional Heterogeneous Data Analysis
# with Generalized Estimating Equations`.


# Before you run the script, make sure you have the geeVerse package installed.
# It can be installed from a source package provided or directly from CRAN.
library(geeVerse)
library(tictoc)
tic()
# Section 4.1 -------------------------------------------------------------

## Example 1 Demonstration Script-----------------------------------------
set.seed(2024)
nsub = 100
sim_data = generateData(nsub= nsub, nobs = rep(10, nsub),  p = 200,
       beta0 = c(rep(1,7),rep(0,193)), rho = 0.6, correlation = "exchangeable")
y = sim_data$y
x = sim_data$X

QPGEE_results = qpgee(x, y, nobs=rep(10, nsub), correlation = "exchangeable",
                      tau = 0.9, max_it=100, method = "HBIC")
names(QPGEE_results)
compile_result(QPGEE_results, beta0= c(rep(1,7),rep(0,193)), cutoff = 0.1)

## Monte Carlo Simulation Script for Table 2 -----------------------------
n_sim = 1 # number of simulation replications
res_list = NULL # initialize a list to store all the fitted models
# start simulation
for (m in 1:n_sim){
  #generate data
  nsub = 100
  set.seed(2024)
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
compile_result(res_list[["independence"]][[as.character(0.9)]],
               beta0= c(rep(1,7),rep(0,193)))


## Example 2 Demonstration Script-----------------------------------------
qpgee(x, y , nobs = rep(1,1000), correlation = "independence", tau = 0.9,
      max_it = 100, method = "HBIC")

## Monte Carlo Simulation  Script for Table 3 -----------------------------
n_sim = 1 # number of simulation replications
res_list = NULL # initialize a list to store all the fitted models
# start simulation
for (m in 1:n_sim){
  #generate data
  nsub = 100
  set.seed(2024)
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


## Example 3 Demonstration Script-----------------------------------------
#generate imbalanced data
nobs = c(rep(10,nsub/2),rep(5,nsub/2))
sim_data = generateData(nsub= 100, nobs = nobs, p = 200,
                        beta0 = c(rep(1,7),rep(0,193)), rho = 0.6, correlation = "exchangeable")
y = sim_data$y
x = sim_data$X
# fit qpgee
qpgee(x, y, nobs = c(rep(10,nsub/2),rep(5,nsub/2)),
      correlation = "exchangeable", tau =0.9, max_it=100, method = "HBIC")

## Monte Carlo Simulation Script for Table 4 ---------------------------------
n_sim = 1 # number of simulation replications
res_list = NULL # initialize a list to store all the fitted models
# start simulation
for (m in 1:n_sim){
  #generate data
  nsub = 100
  set.seed(2024)
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

## Demonstration Script for example 4 ---------------------------
set.seed(2024)
sim_data = generateData(nsub = nsub, nobs = rep(10, nsub), p = 1000,
                        beta0 = c(rep(1,7),rep(0,993)), rho = 0.6, correlation = "exchangeable")
y = sim_data$y
x = sim_data$X
x_sis_ind = SIS::SIS(x,y,nsis=200)$sis.ix0
xsis = x[,x_sis_ind]
QPGEE_results_highdim = qpgee(xsis, y, nobs=rep(10, nsub),
                              correlation = "exchangeable", tau=0.9, max_it=100, method = "HBIC", ncore = 10)

QPGEE_results_highdim$beta <- replace(vector("numeric", 1000), x_sis_ind, QPGEE_results_highdim$beta)
compile_result(QPGEE_results_highdim,beta0= c(rep(1,7),rep(0,993)))


# Section 4.3 -------------------------------------------------------------
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


# Section 5.1 -------------------------------------------------------------

# load data
data("simuGene")
set.seed(2024)
sim_data <- generateData(nsub = 1000, nobs = rep(5,1000), p = 50,
                         beta0 = c(rep(1,9),rep(0,41)), rho = 0.6, correlation = "exchangeable",
                         ka = 0.5, SNPs = simuGene[,1:25])
y = sim_data$y
x = sim_data$X

PQGEE_results_median <- qpgee(x, y, nobs=rep(5,1000), correlation="exchangeable",
                              tau=0.5, intercept= FALSE, method="HBIC", cutoff=10^-4,ncore = 10)
compile_result(PQGEE_results_median,
               beta0= c(rep(1,9),rep(0,41)))



# Section 5.2 -------------------------------------------------------------

## Gene Expression Example: Yeast Cell G1 ---------------------------------
# Load data, the yeastG1 data is provided in geeVerse package
data("yeastG1")
n_obs=4
n_sub=283
x=as.matrix (yeastG1[,c(3:99)])
y=as.matrix(yeastG1$y)
predictors=c("intercept",names(yeastG1[,c(3:99)]))

# obtain beta initial from a cross-sectional non-penalized quantile model.
library(quantreg)
betaint = coefficients(rq( y ~ x, tau = 0.9))
# qpgee with HBIC tuning
PQGEE_results <- qpgee(x, y, nobs = rep(n_obs,n_sub), betaint = betaint,
                       correlation= "unstructured", tau=0.9, intercept=TRUE, method="HBIC", cutoff=10^-4)
# report selected variables
predictors[abs(PQGEE_results$beta)>10^-3]

# PGEE with CV selected penalty level
CVfit_result <- CVfit("y ~ . -id ", id = id, data = yeastG1,
                      lambda.vec = exp(seq(log(10),log(0.1),length.out = 30)), fold = 5)
myfit1 <- PGEE("y ~ . - id ", id = id, data = yeastG1,
               lambda = CVfit_result$lam.opt)
# report selected variables
pgee_index <- which(abs(myfit1$coefficient) > 10^-3)
predictors[pgee_index]


# Section 5.3 -------------------------------------------------------------
# As CCI-779 data is not publicly available, we are not able to provide the script here.
toc()
