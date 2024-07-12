library(geeVerse)
library(quantreg)
data("yeastG1")
nobs=4
nsub<-nrow(yeastG1)/nobs


predictors=c("intercept",names(yeastG1[,c(3:99)]))
x=as.matrix (yeastG1[,c(3:99)])
y=as.matrix(yeastG1$y)

#Apply proposed method with hbic tuning with tau = 0.9
betaint=coefficients(rq(y~0+cbind(1,x),tau = 0.9))
PQGEE_results_high <-  qpgee(x, y, nobs = rep(nobs,nsub),
                               correlation  = "Ind", tau=0.9 , intercept = TRUE, method = "HBIC", cutoff = 10^-4,ncore=10)
predictors[abs(PQGEE_results_high$beta)>10^-3]


#Apply proposed method with hbic tuning with tau = 0.5
betaint=coefficients(rq(y~x,tau = 0.5))
PQGEE_results_median <-  qpgee(x, y, nobs = rep(nobs,nsub), betaint = betaint,
                               correlation  = "Ind", tau=0.5 , intercept = TRUE, method = "HBIC", cutoff = 10^-4, ncore = 10)
predictors[abs(PQGEE_results_median$beta)>10^-3]

#PGEE mean
CVfit_result <- CVfit("y ~ . -id ", id = id, data = yeastG1, lambda.vec = exp(seq(log(10),log(0.1),length.out = 30)), fold = 5)
myfit1 <- PGEE("y ~ . - id ", id = id, data = yeastG1,lambda = CVfit_result$lam.opt)
index1 <- which(abs(myfit1$coefficient) > 10^-3)
predictors[index1]



