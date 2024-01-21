
qpgee_tvec = c()
pgee_own_tvec = c()
pgee_pkg_tvec = c()

for (i in 1:100) {
  #data generation settings
  n=n_sub=200
  p=200
  beta0=rep(1,7)
  p0=length(beta0)
  beta = c(beta0,rep(0,p-p0))
  n_obs<-rep(10,n_sub);
  ka = 1
  rho=0.6
  type="ar"
  dis="normal"

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
  x=X=matrix(rnorm(N*p),N,p)
  y=X%*%beta+(1+ka*abs(X[,1]))*e

  #fit qpgee
  qpgee_time = system.time(qpgee_fit <- qpgee(x,y,tau=0.5,nk=n_obs))[3]
  qpgee_tvec = c(qpgee_tvec,qpgee_time)

  #fit own pgee
  data = data.frame(X,y,id)
  #PGEE_fit = PGEE("y ~.-id-1",id = id, data = data,corstr = "exchangeable",lambda=0.01)
  pgee_own_time = system.time(PGEE_own_fit <- PGEE_own("y ~.-id-1",id = id,data = data,
                                       corstr = "exchangeable",lambda=0.1))[3]
  pgee_own_tvec = c(pgee_own_tvec,pgee_own_time)
  #fit pgee from PGEE package
  pgee_pkg_time = system.time(PGEE_pkg <- PGEE::PGEE("y ~.-id-1",id = id,data = data,
                                     corstr = "exchangeable",lambda=0.1))[3]
  pgee_pkg_tvec = c(pgee_pkg_tvec,pgee_pkg_time)

}

report_sum <- function(x){
  c(mean(x),min(x),max(x),median(x))
}

report_sum(qpgee_tvec)
report_sum(pgee_own_tvec)
report_sum(pgee_pkg_tvec)
qpgee_tvec_p200 = qpgee_tvec
pgee_own_tvec_p200 = pgee_own_tvec
pgee_pkg_tvec_p200 = pgee_pkg_tvec
report_sum(qpgee_tvec_p200)
  report_sum(pgee_own_tvec_p200)
    report_sum(pgee_pkg_tvec_p200)
qpgee_tvec = c()
pgee_own_tvec = c()
pgee_pkg_tvec = c()

for (i in 1:50) {
  #data generation settings
  n=n_sub=200
  p=1000
  beta0=rep(1,7)
  p0=length(beta0)
  beta = c(beta0,rep(0,p-p0))
  n_obs<-rep(10,n_sub);
  ka = 1
  rho=0.6
  type="ar"
  dis="normal"

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
  x=X=matrix(rnorm(N*p),N,p)
  y=X%*%beta+(1+ka*abs(X[,1]))*e

  #fit qpgee
  qpgee_time = system.time(qpgee_fit <- qpgee(x,y,tau=0.5,nk=n_obs))[3]
  qpgee_tvec = c(qpgee_tvec,qpgee_time)

  #fit own pgee
  data = data.frame(X,y,id)
  #PGEE_fit = PGEE("y ~.-id-1",id = id, data = data,corstr = "exchangeable",lambda=0.01)
  pgee_own_time = system.time(PGEE_own_fit <- PGEE_own("y ~.-id-1",id = id,data = data,
                                                       corstr = "exchangeable",lambda=0.1))[3]
  pgee_own_tvec = c(pgee_own_tvec,pgee_own_time)
  #fit pgee from PGEE package
  pgee_pkg_time = system.time(PGEE_pkg <- PGEE::PGEE("y ~.-id-1",id = id,data = data,
                                                     corstr = "exchangeable",lambda=0.1))[3]
  pgee_pkg_tvec = c(pgee_pkg_tvec,pgee_pkg_time)

}

report_sum <- function(x){
  c(mean(x),min(x),max(x),median(x))
}

report_sum(qpgee_tvec)
report_sum(pgee_own_tvec)
report_sum(pgee_pkg_tvec)

