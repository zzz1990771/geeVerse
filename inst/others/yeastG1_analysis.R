library(quantreg)

data("yeastG1")
n_obs=4
n_sub<-nrow(yeastG1)/n_obs


predictors=c("intercept",names(yeastG1[,c(3:99)]))
x=as.matrix (yeastG1[,c(3:99)])
y=as.matrix(yeastG1$y)

betaint=coefficients(rq(y~0+x,tau = 0.9))

#Apply proposed method with hbic tuning with tau = 0.9
PQGEE_results_high <-  qpgee_tune(cbind(1,x), y, nk = rep(n_obs,n_sub), correlation = "Ind",
                                  tau=0.9 ,method = "HBIC",cutoff = 10^-4,ncore = 10)
predictors[PQGEE_results_high$X_selected]


#Apply proposed method with hbic tuning with tau = 0.5
PQGEE_results_median <-  qpgee_tune(cbind(1,x), y, nk = rep(n_obs,n_sub), correlation = "Ind",
                             tau=0.5 ,method = "HBIC",cutoff = 10^-4,ncore = 10)
predictors[PQGEE_results_median$X_selected]


CVfit_result <- CVfit("y ~ . -id ", id = id, data = yeastG1, lambda.vec = exp(seq(log(10),log(0.1),length.out = 30)), fold = 5)
myfit1 <- PGEE("y ~ . - id ", id = id, data = yeastG1,lambda = CVfit_result$lam.opt)
index1 <- which(abs(myfit1$coefficient) > 10^-3)
predictors[index1]



