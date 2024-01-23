####PQIF####
###############################################################################
###############################################################################
## ##
## This is R coding for the paper ##
## "Model selection for correlated data with diverging number of parameters" ##
## by Hyunkeun Cho and Annie Qu ##
## Department of Statistics, University of Illinois at Urbana-Champaign ##
## ##
## For the model selection of the basis matrices for the correlation matrix ##
## (Zhou and Qu, 2012, JASA), see the website for the paper and R coding: ##
## https://netfiles.uiuc.edu/anniequ/public_html/ ##
###############################################################################
###############################################################################
# library(MASS)
# #install.packages('plus')
# library(plus)
# library(geepack)
# library(mvtnorm)

#####################################
#Function for generating binary data#
#####################################

generate_binaryOriginal<-function(n,pn,m,str.y,beta.true,rho,sigma2){
  #n=100;pn=100;m=n_obs;str.y="exch";beta.true=true.alpha;rho=rho
  #n : sample size
  #pn : the number of parameter
  #m : cluster size
  #rho : correlation coefficient
  #str.y : correlation structure
  #beta.true : true parameter
  #sigma2 : null
  x=NULL
  error=NULL
  if(str.y=="ar1"){
    mm <- diag(m)
    for(i in 1:(m-1))
    {
      m1 <- matrix(rep(0, m * m), m)
      for (k in 1:m) {
        for (l in 1:m) {
          if (abs(k - l) == i)
            m1[k, l] <- rho^i
        }
      }
      mm<-mm+m1
    }
  }
  if(str.y=="exch"){
    mm <- diag(m)
    m1 <- matrix(rep(1, m * m), m) - mm
    mm <- mm + rho*m1
  }
  a <- runif(n*m*1,0,0.8)
  a <- matrix(a,nrow=(n*m))
  b <- runif(n*m*2,0,0.8)
  b <- matrix(b,nrow=(n*m))
  c <- runif(n*m*(pn-3),0,1)
  c <- matrix(c,nrow=(n*m))
  x <- as.matrix(cbind(a,b,c))
  u <- x%*%beta.true
  muu <- exp(u)/(1+exp(u))
  mu <- NULL
  for(i in 1:n){
    ep0=ep(mu=c(muu[(m*(i-1)+1):(m*i),]), R=mm, nRep=1) #,muu[ni*(i-1)+2],muu[ni*(i-1)+3],muu[
    muy<-matrix(ep0$y,ncol=1)
    mu<-rbind(mu,muy)
  }
  id.vect<-rep(1:n, each=m)
  Y<-as.vector(mu)
  X<-as.matrix(x)
  return(list(y.vect = Y, x.mat = X, id.vect=id.vect, m=m, n=n, pn=pn, beta.true=beta.true, rho=rho))
}


################
#Logit function#
################
mm <- function(x){
  value <- exp(x) / (1+exp(x))
  value[value == 1] <- 1 - 10^(-6)
  value[value < 1e-156] <- 10^(-6)
  return(value)
}


##############################################
#Estimate the initial estimator of parameters#
##############################################
qifforspline<-function(X,Y,cfix,nsub,ni,working){#a,str.y,
  # X=bs;cfix=Z%*%b;nsub=n_sub;ni=n_obs;working=corrQIF[corrStructure]
  #X=Covbs;nsub=n_sub;ni=n_obs;working=corrQIF[corrStructure]
  #X : explanatory variable
  #Y : response variable
  #nsub : samle size
  #ni : cluster size
  #str.y : correlation structure
  #a : correlation coefficient
  #working : working correlation structure
  X<-as.matrix(X)
  Y<-as.matrix(Y)
  id<-matrix(rep((1:nsub),each=ni),nsub*ni)
  obs <- lapply(split(id, id), "length")
  nobs <- as.numeric(obs) #number of cluster for each subject
  np<-dim(X)[[2]] #number of parameters in the selected marginal model
  #if(fam="gaussian"){
  #G <- geeglm(Y~X-1,id=id)
  #}
  #if(fam="logit"){
  G <- geeglm(Y~X-1, family=binomial(link="logit"),id=id)
  #}


  beta<-as.numeric(G$coefficients)
  betadiff <- 1
  iteration <- 0
  betanew <- beta
  tol<-1e-12
  maxiter<-1000
  if(working=="indep"){
    while (betadiff > tol && iteration < maxiter) {
      beta <- betanew
      sumg <- matrix(rep(0, np), nrow = np)
      sumc <- matrix(rep(0, np * np), nrow = np)
      arsumg <- matrix(rep(0, np), nrow = np)
      arsumc <- matrix(rep(0, np * np), nrow = np)
      gi <- matrix(rep(0, np), nrow = np)
      arsumgfirstdev <- matrix(rep(0, np * np), nrow = np)
      firstdev <- matrix(rep(0, np * np), nrow = np)
      loc1 <- 0
      loc2 <- 0
      for (i in 1:nsub) {
        loc1 <- loc2 + 1
        loc2 <- loc1 + nobs[i] - 1
        yi <- as.matrix(Y[loc1:loc2, ])
        xi <- X[loc1:loc2, ]
        cfixi<-cfix[loc1:loc2]
        ni <- nrow(yi)
        m0 <- diag(ni)
        ui <- mm(xi %*% beta+cfixi)
        fui <- log(ui) - log(1 - ui)
        if(dim(ui)[[1]]==1)
        {
          fui_dev <- diag(ui) %*% diag(1 - ui)
          vui <- diag(sqrt(1/ui)) %*% diag(sqrt(1/(1 - ui)))
          wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
        } #fui_dev %*% vui %*% m0 %*% vui is identity maatrix
        if(dim(ui)[[1]]!=1){
          fui_dev <- diag((as.vector(ui))) %*% diag(as.vector(1 - ui))
          vui <- diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1 - ui))))
          wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
        }
        gi <- (1/nsub) * wi %*% (yi - ui)
        arsumc <- arsumc + gi %*% t(gi)
        arsumg <- arsumg + gi
        di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
        firstdev <- di0
        arsumgfirstdev <- arsumgfirstdev + firstdev
      }
      arcinv = geninv(arsumc)
      Q <- t(arsumg) %*% arcinv %*% arsumg
      arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
      arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev
      invarqif2dev <- geninv(arqif2dev)
      betanew <- beta - invarqif2dev %*% arqif1dev
      betadiff <- abs(sum(betanew - beta))
      iteration <- iteration + 1
    }
  }
  if(working!="indep"){
    while (betadiff > tol && iteration < maxiter) {
      beta <- betanew
      sumg <- matrix(rep(0, 2 * np), nrow = 2 * np)
      sumc <- matrix(rep(0, 2 * np * 2 * np), nrow = 2 * np)
      arsumg <- matrix(rep(0, 2 * np), nrow = 2 * np)
      arsumc <- matrix(rep(0, 2 * np * 2 * np), nrow = 2 * np)
      gi <- matrix(rep(0, 2 * np), nrow = 2 * np)
      arsumgfirstdev <- matrix(rep(0, 2 * np * np), nrow = 2 * np)
      firstdev <- matrix(rep(0, 2 * np * np), nrow = 2 * np)
      loc1 <- 0
      loc2 <- 0
      for (i in 1:nsub) {
        loc1 <- loc2 + 1
        loc2 <- loc1 + nobs[i] - 1
        yi <- as.matrix(Y[loc1:loc2, ])
        xi <- X[loc1:loc2, ]
        cfixi<-cfix[loc1:loc2]
        ni <- nrow(yi)
        m0 <- diag(ni)
        if(working=="ar1"){
          m1 <- matrix(rep(0, ni * ni), ni)
          for (k in 1:ni) {
            for (l in 1:ni) {
              if (abs(k - l) == 1)
                m1[k, l] <- 1
            }
          }
        }
        if(working=="exch"){
          m1 <- matrix(rep(1, ni * ni), ni) - m0
        }
        ui <- mm(xi %*% beta+cfixi)
        fui <- log(ui)-log(1-ui)
        if(dim(ui)[[1]]==1)
        {
          fui_dev <- diag(ui) %*% diag(1 - ui)
          vui <- diag(sqrt(1/ui)) %*% diag(sqrt(1/(1 - ui)))
          wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
          zi <- xi %*% fui_dev %*% vui %*% m1 %*% vui
        }
        if(dim(ui)[[1]]!=1){
          fui_dev <- diag((as.vector(ui))) %*% diag(as.vector(1 - ui))
          vui <- diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1 - ui))))
          wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
          zi <- t(xi) %*% fui_dev %*% vui %*% m1 %*% vui
        }
        gi0 <- (1/nsub) * wi %*% (yi - ui)
        gi1 <- (1/nsub) * zi %*% (yi - ui)
        gi[1:np, ] <- gi0
        gi[(np + 1):(2 * np), ] <- gi1
        arsumc <- arsumc + gi %*% t(gi)
        arsumg <- arsumg + gi
        di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
        di1 <- -(1/nsub) * zi %*% fui_dev %*% xi
        firstdev[1:np, ] <- di0
        firstdev[(np + 1):(2 * np), ] <- di1
        arsumgfirstdev <- arsumgfirstdev + firstdev
      }
      arcinv = geninv(arsumc)
      Q <- t(arsumg) %*% arcinv %*% arsumg
      arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
      arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev
      invarqif2dev <- geninv(arqif2dev)
      betanew <- beta - invarqif2dev %*% arqif1dev
      betadiff <- t(betanew - beta)%*%(betanew - beta)
      iteration <- iteration + 1
      #print("betanew from iqif function")
      #print(betanew)
    }
  }
  return(list(betanew=betanew))

}
iqif<-function(X,Y,nsub,ni,working){#a,str.y,
  #X=Covbs;nsub=n_sub;ni=n_obs;working=corrQIF[corrStructure]
  #X : explanatory variable
  #Y : response variable
  #nsub : samle size
  #ni : cluster size
  #str.y : correlation structure
  #a : correlation coefficient
  #working : working correlation structure
  X<-as.matrix(X)
  Y<-as.matrix(Y)
  id<-matrix(rep((1:nsub),each=ni),nsub*ni)
  obs <- lapply(split(id, id), "length")
  nobs <- as.numeric(obs) #number of cluster for each subject
  np<-dim(X)[[2]] #number of parameters in the selected marginal model
  #if(fam="gaussian"){
  #G <- geeglm(Y~X-1,id=id)
  #}
  #if(fam="logit"){
  G <- geeglm(Y~X-1, family=binomial(link="logit"),id=id)
  #}


  beta<-as.numeric(G$coefficients)
  betadiff <- 1
  iteration <- 0
  betanew <- beta
  tol<-1e-12
  maxiter<-1000
  if(working=="indep"){
    while (betadiff > tol && iteration < maxiter) {
      beta <- betanew
      sumg <- matrix(rep(0, np), nrow = np)
      sumc <- matrix(rep(0, np * np), nrow = np)
      arsumg <- matrix(rep(0, np), nrow = np)
      arsumc <- matrix(rep(0, np * np), nrow = np)
      gi <- matrix(rep(0, np), nrow = np)
      arsumgfirstdev <- matrix(rep(0, np * np), nrow = np)
      firstdev <- matrix(rep(0, np * np), nrow = np)
      loc1 <- 0
      loc2 <- 0
      for (i in 1:nsub) {
        loc1 <- loc2 + 1
        loc2 <- loc1 + nobs[i] - 1
        yi <- as.matrix(Y[loc1:loc2, ])
        xi <- X[loc1:loc2, ]
        ni <- nrow(yi)
        m0 <- diag(ni)
        ui <- mm(xi %*% beta)
        fui <- log(ui) - log(1 - ui)
        if(dim(ui)[[1]]==1)
        {
          fui_dev <- diag(ui) %*% diag(1 - ui)
          vui <- diag(sqrt(1/ui)) %*% diag(sqrt(1/(1 - ui)))
          wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
        } #fui_dev %*% vui %*% m0 %*% vui is identity maatrix
        if(dim(ui)[[1]]!=1){
          fui_dev <- diag((as.vector(ui))) %*% diag(as.vector(1 - ui))
          vui <- diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1 - ui))))
          wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
        }
        gi <- (1/nsub) * wi %*% (yi - ui)
        arsumc <- arsumc + gi %*% t(gi)
        arsumg <- arsumg + gi
        di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
        firstdev <- di0
        arsumgfirstdev <- arsumgfirstdev + firstdev
      }
      arcinv = geninv(arsumc)
      Q <- t(arsumg) %*% arcinv %*% arsumg
      arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
      arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev
      invarqif2dev <- geninv(arqif2dev)
      betanew <- beta - invarqif2dev %*% arqif1dev
      betadiff <- abs(sum(betanew - beta))
      iteration <- iteration + 1
    }
  }
  if(working!="indep"){
    while (betadiff > tol && iteration < maxiter) {
      beta <- betanew
      sumg <- matrix(rep(0, 2 * np), nrow = 2 * np)
      sumc <- matrix(rep(0, 2 * np * 2 * np), nrow = 2 * np)
      arsumg <- matrix(rep(0, 2 * np), nrow = 2 * np)
      arsumc <- matrix(rep(0, 2 * np * 2 * np), nrow = 2 * np)
      gi <- matrix(rep(0, 2 * np), nrow = 2 * np)
      arsumgfirstdev <- matrix(rep(0, 2 * np * np), nrow = 2 * np)
      firstdev <- matrix(rep(0, 2 * np * np), nrow = 2 * np)
      loc1 <- 0
      loc2 <- 0
      for (i in 1:nsub) {
        loc1 <- loc2 + 1
        loc2 <- loc1 + nobs[i] - 1
        yi <- as.matrix(Y[loc1:loc2, ])
        xi <- X[loc1:loc2, ]
        ni <- nrow(yi)
        m0 <- diag(ni)
        if(working=="ar1"){
          m1 <- matrix(rep(0, ni * ni), ni)
          for (k in 1:ni) {
            for (l in 1:ni) {
              if (abs(k - l) == 1)
                m1[k, l] <- 1
            }
          }
        }
        if(working=="exch"){
          m1 <- matrix(rep(1, ni * ni), ni) - m0
        }
        ui <- mm(xi %*% beta)
        fui <- log(ui)-log(1-ui)
        if(dim(ui)[[1]]==1)
        {
          fui_dev <- diag(ui) %*% diag(1 - ui)
          vui <- diag(sqrt(1/ui)) %*% diag(sqrt(1/(1 - ui)))
          wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
          zi <- xi %*% fui_dev %*% vui %*% m1 %*% vui
        }
        if(dim(ui)[[1]]!=1){
          fui_dev <- diag((as.vector(ui))) %*% diag(as.vector(1 - ui))
          vui <- diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1 - ui))))
          wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
          zi <- t(xi) %*% fui_dev %*% vui %*% m1 %*% vui
        }
        gi0 <- (1/nsub) * wi %*% (yi - ui)
        gi1 <- (1/nsub) * zi %*% (yi - ui)
        gi[1:np, ] <- gi0
        gi[(np + 1):(2 * np), ] <- gi1
        arsumc <- arsumc + gi %*% t(gi)
        arsumg <- arsumg + gi
        di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
        di1 <- -(1/nsub) * zi %*% fui_dev %*% xi
        firstdev[1:np, ] <- di0
        firstdev[(np + 1):(2 * np), ] <- di1
        arsumgfirstdev <- arsumgfirstdev + firstdev
      }
      arcinv = geninv(arsumc)
      Q <- t(arsumg) %*% arcinv %*% arsumg
      arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
      arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev
      invarqif2dev <- geninv(arqif2dev)
      betanew <- beta - invarqif2dev %*% arqif1dev
      betadiff <- t(betanew - beta)%*%(betanew - beta)
      iteration <- iteration + 1
      #print("betanew from iqif function")
      #print(betanew)
    }
  }
  return(list(betanew=betanew))

}

igausqif<-function(X,Y,nsub,ni,working){#a,str.y,

  #X=cbind(X0,Z0);nsub=n_sub;ni=n_obs; working=corrQIF[corrStructure]
  #X : explanatory variable
  #Y : response variable
  #nsub : samle size
  #ni : cluster size
  #str.y : correlation structure
  #a : correlation coefficient
  #working : working correlation structure
  X<-as.matrix(X)
  Y<-as.matrix(Y)
  id<-matrix(rep((1:nsub),each=ni),nsub*ni)
  obs <- lapply(split(id, id), "length")
  nobs <- as.numeric(obs) #number of cluster for each subject
  np<-dim(X)[[2]] #number of parameters in the selected marginal model
  #if(fam="gaussian"){
  #G <- geeglm(Y~X-1,id=id)
  #}
  #if(fam="logit"){
  G <- geeglm(Y~X-1,id=id)
  #}


  beta<-as.numeric(G$coefficients)
  betadiff <- 1
  iteration <- 0
  betanew <- beta
  tol<-1e-12
  maxiter<-1000
  if(working=="indep"){
    while (betadiff > tol && iteration < maxiter) {
      beta <- betanew
      sumg <- matrix(rep(0, np), nrow = np)
      sumc <- matrix(rep(0, np * np), nrow = np)
      arsumg <- matrix(rep(0, np), nrow = np)
      arsumc <- matrix(rep(0, np * np), nrow = np)
      gi <- matrix(rep(0, np), nrow = np)
      arsumgfirstdev <- matrix(rep(0, np * np), nrow = np)
      firstdev <- matrix(rep(0, np * np), nrow = np)
      loc1 <- 0
      loc2 <- 0
      for (i in 1:nsub) {
        loc1 <- loc2 + 1
        loc2 <- loc1 + nobs[i] - 1
        yi <- as.matrix(Y[loc1:loc2, ])
        xi <- X[loc1:loc2, ]
        ni <- nrow(yi)
        m0 <- diag(ni)

        ui <- xi %*% beta
        fui <- ui
        # ui <- mm(xi %*% beta)
        #fui <- log(ui) - log(1 - ui)
        if(dim(ui)[[1]]==1)
        {
          fui_dev <- diag(ni)
          vui <-diag(ni)
          wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
        } #fui_dev %*% vui %*% m0 %*% vui is identity maatrix
        if(dim(ui)[[1]]!=1){
          fui_dev <- diag(ni)#diag((as.vector(ui))) %*% diag(as.vector(1 - ui))
          vui <- diag(ni)#diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1 - ui))))
          wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
        }
        gi <- (1/nsub) * wi %*% (yi - ui)
        arsumc <- arsumc + gi %*% t(gi)
        arsumg <- arsumg + gi
        di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
        firstdev <- di0
        arsumgfirstdev <- arsumgfirstdev + firstdev
      }
      arcinv = geninv(arsumc)
      Q <- t(arsumg) %*% arcinv %*% arsumg
      arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
      arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev
      invarqif2dev <- geninv(arqif2dev)
      betanew <- beta - invarqif2dev %*% arqif1dev
      betadiff <- abs(sum(betanew - beta))
      iteration <- iteration + 1
    }
  }
  if(working!="indep"){
    while (betadiff > tol && iteration < maxiter) {
      beta <- betanew
      sumg <- matrix(rep(0, 2 * np), nrow = 2 * np)
      sumc <- matrix(rep(0, 2 * np * 2 * np), nrow = 2 * np)
      arsumg <- matrix(rep(0, 2 * np), nrow = 2 * np)
      arsumc <- matrix(rep(0, 2 * np * 2 * np), nrow = 2 * np)
      gi <- matrix(rep(0, 2 * np), nrow = 2 * np)
      arsumgfirstdev <- matrix(rep(0, 2 * np * np), nrow = 2 * np)
      firstdev <- matrix(rep(0, 2 * np * np), nrow = 2 * np)
      loc1 <- 0
      loc2 <- 0
      for (i in 1:nsub) {
        loc1 <- loc2 + 1
        loc2 <- loc1 + nobs[i] - 1
        yi <- as.matrix(Y[loc1:loc2, ])
        xi <- X[loc1:loc2, ]
        ni <- nrow(yi)
        m0 <- diag(ni)
        if(working=="ar1"){
          m1 <- matrix(rep(0, ni * ni), ni)
          for (k in 1:ni) {
            for (l in 1:ni) {
              if (abs(k - l) == 1)
                m1[k, l] <- 1
            }
          }
        }
        if(working=="exch"){
          m1 <- matrix(rep(1, ni * ni), ni) - m0
        }
        ui <- xi %*% beta
        fui <- ui
        if(dim(ui)[[1]]==1)
        {
          fui_dev <- diag(ni)
          vui <- diag(ni)
          wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
          zi <- xi %*% fui_dev %*% vui %*% m1 %*% vui
        }
        if(dim(ui)[[1]]!=1){
          fui_dev <- diag(ni)
          vui <- diag(ni)
          wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
          zi <- t(xi) %*% fui_dev %*% vui %*% m1 %*% vui
        }
        gi0 <- (1/nsub) * wi %*% (yi - ui)
        gi1 <- (1/nsub) * zi %*% (yi - ui)
        gi[1:np, ] <- gi0
        gi[(np + 1):(2 * np), ] <- gi1
        arsumc <- arsumc + gi %*% t(gi)
        arsumg <- arsumg + gi
        di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
        di1 <- -(1/nsub) * zi %*% fui_dev %*% xi
        firstdev[1:np, ] <- di0
        firstdev[(np + 1):(2 * np), ] <- di1
        arsumgfirstdev <- arsumgfirstdev + firstdev
      }
      arcinv = geninv(arsumc)
      Q <- t(arsumg) %*% arcinv %*% arsumg
      arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
      arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev
      invarqif2dev <- geninv(arqif2dev)
      betanew <- beta - invarqif2dev %*% arqif1dev
      betadiff <- t(betanew - beta)%*%(betanew - beta)
      iteration <- iteration + 1
    }
  }
  return(list(betanew=betanew))

}


########################################
#Estimate the parameter of oracle model#
########################################
oqif<-function(nsub,ni,X,Y,working){
  X<-as.matrix(X)
  Y<-as.matrix(Y)
  id<-matrix(rep((1:nsub),each=ni),nsub*ni)
  obs <- lapply(split(id, id), "length")
  nobs <- as.numeric(obs) #number of cluster for each subject
  np<-dim(X)[[2]] #number of parameters in the selected marginal model
  G <- geeglm(Y~X-1, family=binomial(link="logit"),id=id)
  beta<-as.numeric(G$coefficients)
  betadiff <- 1
  iteration <- 0
  betanew <- beta
  tol<-1e-12
  maxiter<-1000
  if(working=="indep"){while (betadiff > tol && iteration < maxiter) {
    beta <- betanew
    sumg <- matrix(rep(0, np), nrow = np)
    sumc <- matrix(rep(0, np * np), nrow = np)
    arsumg <- matrix(rep(0, np), nrow = np)
    arsumc <- matrix(rep(0, np * np), nrow = np)
    gi <- matrix(rep(0, np), nrow = np)
    arsumgfirstdev <- matrix(rep(0, np * np), nrow = np)
    firstdev <- matrix(rep(0, np * np), nrow = np)
    loc1 <- 0
    loc2 <- 0
    for (i in 1:nsub) {
      loc1 <- loc2 + 1
      loc2 <- loc1 + nobs[i] - 1
      yi <- as.matrix(Y[loc1:loc2, ])
      xi <- X[loc1:loc2, ]
      ni <- nrow(yi)
      m0 <- diag(ni)
      ui <- mm(xi %*% beta)
      fui <- log(ui) - log(1 - ui)
      if(dim(ui)[[1]]==1)
      {
        fui_dev <- diag(ui) %*% diag(1 - ui)
        vui <- diag(sqrt(1/ui)) %*% diag(sqrt(1/(1 - ui)))
        wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
      } #fui_dev %*% vui %*% m0 %*% vui is identity maatrix
      if(dim(ui)[[1]]!=1){
        fui_dev <- diag((as.vector(ui))) %*% diag(as.vector(1 - ui))
        vui <- diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1 - ui))))
        wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
      }
      gi <- (1/nsub) * wi %*% (yi - ui)
      arsumc <- arsumc + gi %*% t(gi)
      arsumg <- arsumg + gi
      di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
      firstdev <- di0
      arsumgfirstdev <- arsumgfirstdev + firstdev
    }
    arcinv = geninv(arsumc)
    Q <- t(arsumg) %*% arcinv %*% arsumg
    arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
    arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev
    invarqif2dev <- geninv(arqif2dev)
    betanew <- beta - invarqif2dev %*% arqif1dev
    betadiff <- abs(sum(betanew - beta))
    iteration <- iteration + 1
  }
  }
  if(working!="indep"){
    while (betadiff > tol && iteration < maxiter) {
      beta <- betanew
      sumg <- matrix(rep(0, 2 * np), nrow = 2 * np)
      sumc <- matrix(rep(0, 2 * np * 2 * np), nrow = 2 * np)
      arsumg <- matrix(rep(0, 2 * np), nrow = 2 * np)
      arsumc <- matrix(rep(0, 2 * np * 2 * np), nrow = 2 * np)
      gi <- matrix(rep(0, 2 * np), nrow = 2 * np)
      arsumgfirstdev <- matrix(rep(0, 2 * np * np), nrow = 2 * np)
      firstdev <- matrix(rep(0, 2 * np * np), nrow = 2 * np)
      loc1 <- 0
      loc2 <- 0
      for (i in 1:nsub) {
        loc1 <- loc2 + 1
        loc2 <- loc1 + nobs[i] - 1
        yi <- as.matrix(Y[loc1:loc2, ])
        xi <- X[loc1:loc2, ]
        ni <- nrow(yi)
        m0 <- diag(ni)
        if(working=="ar1"){
          m1 <- matrix(rep(0, ni * ni), ni)
          for (k in 1:ni) {
            for (l in 1:ni) {
              if (abs(k - l) == 1)
                m1[k, l] <- 1
            }
          }
        }
        if(working=="exch"){
          m1 <- matrix(rep(1, ni * ni), ni) - m0
        }
        ui <- mm(xi %*% beta)
        fui <- log(ui)-log(1-ui)
        if(dim(ui)[[1]]==1)
        {
          fui_dev <- diag(ui) %*% diag(1 - ui)
          vui <- diag(sqrt(1/ui)) %*% diag(sqrt(1/(1 - ui)))
          wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
          zi <- xi %*% fui_dev %*% vui %*% m1 %*% vui
        }
        if(dim(ui)[[1]]!=1){
          fui_dev <- diag((as.vector(ui))) %*% diag(as.vector(1 - ui))
          vui <- diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1 - ui))))
          wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
          zi <- t(xi) %*% fui_dev %*% vui %*% m1 %*% vui
        }
        gi0 <- (1/nsub) * wi %*% (yi - ui)
        gi1 <- (1/nsub) * zi %*% (yi - ui)
        gi[1:np, ] <- gi0
        gi[(np + 1):(2 * np), ] <- gi1
        arsumc <- arsumc + gi %*% t(gi)
        arsumg <- arsumg + gi
        di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
        di1 <- -(1/nsub) * zi %*% fui_dev %*% xi
        firstdev[1:np, ] <- di0
        firstdev[(np + 1):(2 * np), ] <- di1
        arsumgfirstdev <- arsumgfirstdev + firstdev
      }
      arcinv = geninv(arsumc)
      Q <- t(arsumg) %*% arcinv %*% arsumg
      arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
      arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev
      invarqif2dev <- geninv(arqif2dev)
      betanew <- beta - invarqif2dev %*% arqif1dev
      betadiff <- t(betanew - beta)%*%(betanew - beta)
      iteration <- iteration + 1
    }
  }
  return(list(betanew=betanew))
}

#######################
#SCAD penalty function#
#######################
pp <- function(z,lambda_si,lambda_pl,a,dsi) {
  #as.numeric(pp(beta,lambda,T)/abs(beta))#scad penalty
  #print(sprintf("SI dimension is %g",dsi))
  x <- abs(z)
  x <- matrix(x, ncol=1)
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
pp_bin <- function(z,lambda_si,lambda_pl,a,dsi) {
  #as.numeric(pp(beta,lambda,T)/abs(beta))#scad penalty
  dsi=dsi #add one because of constant at beginning
  #print(sprintf("SI dimension is %g",dsi))
  x <- abs(z)
  #print("scad coeff")
  # print(z)
  x <- matrix(x, ncol=1)
  value <- matrix(c(rep(0,dim(x)[[1]])),ncol=1)
  #penalize separately here
  for(i in 1:dim(x)[[1]]){
    lambda<-ifelse(i<=dsi,lambda_si,lambda_pl)
    #print(sprintf("coeff is %g and lambda is %g",i,lambda))

    if(x[i] < lambda){value[i] <- lambda}
    else if (x[i] < a*lambda){
      value[i] <- (a*lambda-x[i])/(a-1)}
    else{
      value[i] <- 0}
  }
  return(value)
}
pp_1lambda <- function(z,lambda,a) {
  #as.numeric(pp(beta,lambda,T)/abs(beta))#scad penalty
  x <- abs(z)
  x <- matrix(x, ncol=1)
  value <- matrix(c(rep(0,dim(x)[[1]])),ncol=1)
  #penalize separately here
  for(i in 1:dim(x)[[1]]){
    if(x[i] < lambda){value[i] <- lambda}
    else if (x[i] < a*lambda){
      value[i] <- (a*lambda-x[i])/(a-1)}
    else{
      value[i] <- 0}
  }
  return(value)
}

#######################
#LASSO penalty function#
#######################
pp_Lasso <- function(z,lambda_si,lambda_pl,a,dsi) {
  #as.numeric(pp(beta,lambda,T)/abs(beta))#scad penalty
  #add one because of constant at beginning
  #print(sprintf("SI dimension is %g",dsi))
  x <- abs(z)
  #print("scad coeff")
  # print(z)
  x <- matrix(x, ncol=1)
  value <- matrix(c(rep(0,dim(x)[[1]])),ncol=1)
  #penalize separately here
  for(i in 1:dim(x)[[1]]){
    lambda<-ifelse(i<=dsi,lambda_si,lambda_pl)
    #print(sprintf("coeff is %g and lambda is %g",i,lambda))
    value[i] <- lambda

  }
  return(value)
}
pplasso_1lambda <- function(z,y,a) {
  #as.numeric(pp(beta,lambda,T)/abs(beta))#scad penalty
  x <- abs(z)
  x <- matrix(x, ncol=1)
  value <- matrix(c(rep(0,dim(x)[[1]])),ncol=1)
  for(i in 1:dim(x)[[1]]){
    value[i] <- y
  }
  return(value)
}

######################################################
#Estimate parameter through our proposed method using#
#a local quadratic approximation #
######################################################
sqif<-function(best,X,Y,cfix,lambda_si,lambda_pl,working,dsi,nsub,ni,T){
  # X=Cov;cfix=constant;working=corrQIF[corrStructure];nsub=n_sub;ni=n_obs
  #best: initial estimate
  #lambda : tuning parameter 'lambda'
  #T : tuning parameter 'a'
  #nsub : samle size
  #ni : cluster size
  #X : explanatory variable
  #Y : response variable
  #working : working correlation structure
  X<-as.matrix(X)
  Y<-as.matrix(Y)
  np<-dim(X)[[2]]
  id<-matrix(rep((1:nsub),each=ni),nsub*ni)
  obs <- lapply(split(id, id), "length")
  nobs <- as.numeric(obs) #number of cluster for each subject
  nsub <- length(nobs) #number of subject
  beta<-best#contained the intercept
  betadiff <- 1
  iteration <- 0
  betanew <- beta
  QIF=NULL
  tol<-1e-10#very important this value have to be larger than simulation tol
  maxiter<-1000
  BEST <- as.matrix(rep(1,np),ncol=1)
  QIF=0
  if(working=="indep"){
    while (betadiff > tol && iteration < maxiter && sum(betanew)!=0) {
      beta <- betanew
      sumg <- matrix(rep(0, np), nrow = np)
      sumc <- matrix(rep(0, np * np), nrow = np)
      arsumg <- matrix(rep(0, np), nrow = np)
      arsumc <- matrix(rep(0, np * np), nrow = np)
      gi <- matrix(rep(0, np), nrow = np)
      arsumgfirstdev <- matrix(rep(0, np * np), nrow = np)
      firstdev <- matrix(rep(0, np * np), nrow = np)
      loc1 <- 0
      loc2 <- 0
      for (i in 1:nsub) {
        loc1 <- loc2 + 1
        loc2 <- loc1 + nobs[i] - 1
        yi <- as.matrix(Y[loc1:loc2, ])
        xi <- X[loc1:loc2, ]
        cfixi<-cfix[loc1:loc2]
        ni <- nrow(yi)
        m0 <- diag(ni)
        if(length(beta)!=1){
          ui <- mm(xi %*% beta+cfixi)}#xi %*% beta+cfix
        else{ ui <- mm(t(beta %*% t(xi))+cfixi)}
        fui <- log(ui)-log(1-ui) #link
        if(dim(ui)[[1]]==1)
        {
          fui_dev <- diag(ui) %*% diag(1 - ui)
          vui <- diag(sqrt(1/ui)) %*% diag(sqrt(1/(1 - ui)))
          wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
        } #fui_dev %*% vui %*% m0 %*% vui is identity maatrix
        if(dim(ui)[[1]]!=1){
          fui_dev <- diag((as.vector(ui))) %*% diag(as.vector(1 - ui))
          vui <- diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1 - ui))))
          wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
        }
        gi <- (1/nsub) * wi %*% (yi - ui)
        arsumc <- arsumc + gi %*% t(gi)
        arsumg <- arsumg + gi
        di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
        firstdev <- di0
        arsumgfirstdev <- arsumgfirstdev + firstdev
      } #end of looping through all subjects
      arcinv = geninv(arsumc)
      Q <- t(arsumg) %*% arcinv %*% arsumg
      QIF<-Q
      arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
      arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev
      pe<-as.numeric(pp_bin(beta,lambda_si,lambda_pl,T,dsi)/abs(beta))#scad penalty z,lambda_si,lambda_pl,a,dsi
      if(length(pe)!=1){
        sigma<- nsub*diag(pe)}
      else{sigma <- nsub*pe}
      #print("geninv(arqif2dev + sigma)")
      #print(arqif2dev + sigma)
      invarqif2dev <- geninv(arqif2dev + sigma)
      betanew <- beta - invarqif2dev %*% (arqif1dev + sigma %*% beta)
      a <- NULL
      b <- NULL
      for(i in 1:dim(betanew)[[1]]){
        if (abs(betanew[i])<1e-4) {
          betanew[i] <- 0 }
        else {
          a<-cbind(a,X[,i])
          b<-cbind(b,betanew[i])
        }
      }
      iid <- (1:length(BEST))[BEST!=0]
      BEST[iid] <- betanew
      betadiff <- t(betanew - beta)%*%(betanew - beta)
      if(sum(betanew)!=0){
        X <- as.matrix(a,nrow=nsub*ni)}
      np <- dim(X)[[2]]
      betanew <- as.numeric(b)
      iteration <- iteration + 1
    }
  }
  if(working!="indep"){
    while (betadiff > tol && iteration < maxiter && sum(betanew)!=0) {

      beta <- betanew
      print(beta)
      sumg <- matrix(rep(0, 2 * np), nrow = 2 * np)
      sumc <- matrix(rep(0, 2 * np * 2 * np), nrow = 2 * np)
      arsumg <- matrix(rep(0, 2 * np), nrow = 2 * np)
      arsumc <- matrix(rep(0, 2 * np * 2 * np), nrow = 2 * np)
      gi <- matrix(rep(0, 2 * np), nrow = 2 * np)
      arsumgfirstdev <- matrix(rep(0, 2 * np * np), nrow = 2 * np)
      firstdev <- matrix(rep(0, 2 * np * np), nrow = 2 * np)
      loc1 <- 0
      loc2 <- 0
      for (i in 1:nsub) {
        loc1 <- loc2 + 1
        loc2 <- loc1 + nobs[i] - 1
        yi <- as.matrix(Y[loc1:loc2, ])
        xi <- X[loc1:loc2, ]
        ni <- nrow(yi)
        cfixi<-cfix[loc1:loc2]
        m0 <- diag(ni)
        if(working=="ar1"){
          m1 <- matrix(rep(0, ni * ni), ni)
          for (k in 1:ni) {
            for (l in 1:ni) {
              if (abs(k - l) == 1)
                m1[k, l] <- 1
            }
          }
        }
        if(working=="exch"){
          m1 <- matrix(rep(1, ni * ni), ni) - m0
        }
        if(length(beta)!=1){
          ui <- mm(xi %*% beta+cfixi)
        }else{ ui <- mm(t(beta %*% t(xi))+cfixi)}
        fui <- log(ui)-log(1-ui)
        if(dim(ui)[[1]]==1)
        {
          vui <- diag(sqrt(1/ui)) %*% diag(sqrt(1/(1 - ui)))
          wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
          zi <- xi %*% fui_dev %*% vui %*% m1 %*% vui
        }
        if(dim(ui)[[1]]!=1){
          fui_dev <- diag((as.vector(ui))) %*% diag(as.vector(1 - ui))
          vui <- diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1 - ui))))
          wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
          zi <- t(xi) %*% fui_dev %*% vui %*% m1 %*% vui
        }
        gi0 <- (1/nsub) * wi %*% (yi - ui)
        gi1 <- (1/nsub) * zi %*% (yi - ui)
        # print("wi")
        # print(wi)
        # print("zi")
        # print(zi)
        gi[1:np, ] <- gi0
        gi[(np + 1):(2 * np), ] <- gi1
        arsumc <- arsumc + gi %*% t(gi)
        arsumg <- arsumg + gi
        di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
        di1 <- -(1/nsub) * zi %*% fui_dev %*% xi
        firstdev[1:np, ] <- di0
        firstdev[(np + 1):(2 * np), ] <- di1
        arsumgfirstdev <- arsumgfirstdev + firstdev
      }
      #print("***************gscore")
      #print(dim(arsumg))
      #print("dimbeta")
      # print(length(betanew))
      #  print(dim(arsumc))
      #  arsumc=arsumc+0.0001*diag(nrow(arsumc))
      #print(sprintf("det of weight matrix W %g",det(arsumc)))
      #print("Number of NAs in W inverse is ...")
      #print(sum(is.na(geninv(arsumc))))

      arcinv = geninv(arsumc)

      Q <- t(arsumg) %*% arcinv %*% arsumg #Qif
      #print(Q-QIF)
      QIF<-Q
      arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
      arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev
      pe<-as.numeric(pp_bin(beta,lambda_si,lambda_pl,T,dsi)/abs(beta))#scad penalty
      print(pe)
      if(length(pe)!=1){
        sigma<- nsub*diag(pe)
      }else{sigma <- nsub*pe}
      #print("geninv(arqif2dev + sigma)")
      #print(arqif2dev + sigma)
      invarqif2dev <- geninv(arqif2dev + sigma)
      betanew <- beta - invarqif2dev %*% (arqif1dev + sigma %*% beta)
      print(betanew)
      a <- NULL
      b <- NULL
      for(i in 1:dim(betanew)[[1]]){ #cutoff
        #print("betanew")
        #print(betanew)
        if (abs(betanew[i])<1e-4) {
          betanew[i] <- 0 }
        else {
          a<-cbind(a,X[,i])
          b<-cbind(b,betanew[i])
        }
      }
      iid <- (1:length(BEST))[BEST!=0]
      BEST[iid] <- betanew
      betadiff <- t(betanew - beta)%*%(betanew - beta)
      # print("betadiff is ")
      #  print(betadiff)
      #print("betanew is ")
      #print(betanew)
      # print("QIF")
      #print(Q)
      if(sum(betanew)!=0){
        X <- as.matrix(a,nrow=nsub*ni)}
      np <- dim(X)[[2]]
      betanew <- as.numeric(b)
      iteration <- iteration + 1
      # print(iteration) shut off for now
    }
  }

  return(list(betanew=BEST,QIF=QIF))
}

sqif_LASSO<-function(best,X,Y,cfix,lambda_si,lambda_pl,working,dsi,nsub,ni,T){
  # X=Cov;cfix=constant;working=corrQIF[corrStructure];nsub=n_sub;ni=n_obs
  #best: initial estimate
  #lambda : tuning parameter 'lambda'
  #T : tuning parameter 'a'
  #nsub : samle size
  #ni : cluster size
  #X : explanatory variable
  #Y : response variable
  #working : working correlation structure
  X<-as.matrix(X)
  Y<-as.matrix(Y)
  np<-dim(X)[[2]]
  id<-matrix(rep((1:nsub),each=ni),nsub*ni)
  obs <- lapply(split(id, id), "length")
  nobs <- as.numeric(obs) #number of cluster for each subject
  nsub <- length(nobs) #number of subject
  beta<-best#contained the intercept
  betadiff <- 1
  iteration <- 0
  betanew <- beta
  QIF=NULL
  tol<-1e-10#very important this value have to be larger than simulation tol
  maxiter<-1000
  BEST <- as.matrix(rep(1,np),ncol=1)
  QIF=0
  if(working=="indep"){
    while (betadiff > tol && iteration < maxiter && sum(betanew)!=0) {
      beta <- betanew
      sumg <- matrix(rep(0, np), nrow = np)
      sumc <- matrix(rep(0, np * np), nrow = np)
      arsumg <- matrix(rep(0, np), nrow = np)
      arsumc <- matrix(rep(0, np * np), nrow = np)
      gi <- matrix(rep(0, np), nrow = np)
      arsumgfirstdev <- matrix(rep(0, np * np), nrow = np)
      firstdev <- matrix(rep(0, np * np), nrow = np)
      loc1 <- 0
      loc2 <- 0
      for (i in 1:nsub) {
        loc1 <- loc2 + 1
        loc2 <- loc1 + nobs[i] - 1
        yi <- as.matrix(Y[loc1:loc2, ])
        xi <- X[loc1:loc2, ]
        cfixi<-cfix[loc1:loc2]
        ni <- nrow(yi)
        m0 <- diag(ni)
        # print("cfix")
        # print(cfix)
        if(length(beta)!=1){
          ui <- mm(xi %*% beta+cfixi)}#xi %*% beta+cfix
        else{ ui <- mm(t(beta %*% t(xi))+cfixi)}
        fui <- log(ui)-log(1-ui) #link
        if(dim(ui)[[1]]==1)
        {
          fui_dev <- diag(ui) %*% diag(1 - ui)
          vui <- diag(sqrt(1/ui)) %*% diag(sqrt(1/(1 - ui)))
          wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
        } #fui_dev %*% vui %*% m0 %*% vui is identity maatrix
        if(dim(ui)[[1]]!=1){
          fui_dev <- diag((as.vector(ui))) %*% diag(as.vector(1 - ui))
          vui <- diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1 - ui))))
          wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
        }
        gi <- (1/nsub) * wi %*% (yi - ui)
        arsumc <- arsumc + gi %*% t(gi)
        arsumg <- arsumg + gi
        di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
        firstdev <- di0
        arsumgfirstdev <- arsumgfirstdev + firstdev
      } #end of looping through all subjects
      #print(sprintf("det of weight matrix W %g",det(arsumc)))
      #print("Number of NAs in W inverse is ...")
      # print(sum(is.na(geninv(arsumc))))
      arcinv = geninv(arsumc)
      Q <- t(arsumg) %*% arcinv %*% arsumg
      QIF<-Q
      arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
      arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev
      pe<-as.numeric(pp_Lasso(beta,lambda_si,lambda_pl,T,dsi)/abs(beta))#scad penalty z,lambda_si,lambda_pl,a,dsi
      if(length(pe)!=1){
        sigma<- nsub*diag(pe)}
      else{sigma <- nsub*pe}
      #print("geninv(arqif2dev + sigma)")
      #print(arqif2dev + sigma)
      invarqif2dev <- geninv(arqif2dev + sigma)
      betanew <- beta - invarqif2dev %*% (arqif1dev + sigma %*% beta)
      a <- NULL
      b <- NULL
      for(i in 1:dim(betanew)[[1]]){
        if (abs(betanew[i])<1e-4) {
          betanew[i] <- 0 }
        else {
          a<-cbind(a,X[,i])
          b<-cbind(b,betanew[i])
        }
      }
      iid <- (1:length(BEST))[BEST!=0]
      BEST[iid] <- betanew
      betadiff <- t(betanew - beta)%*%(betanew - beta)
      # print(sprintf("betadiff is %g tol is %g",betadiff,tol))
      # print(betadiff)
      #print("betanew in sqif is ")
      #print(betanew)
      if(sum(betanew)!=0){
        X <- as.matrix(a,nrow=nsub*ni)}
      np <- dim(X)[[2]]
      betanew <- as.numeric(b)
      iteration <- iteration + 1
    }
  }
  if(working!="indep"){
    while (betadiff > tol && iteration < maxiter && sum(betanew)!=0) {
      beta <- betanew
      sumg <- matrix(rep(0, 2 * np), nrow = 2 * np)
      sumc <- matrix(rep(0, 2 * np * 2 * np), nrow = 2 * np)
      arsumg <- matrix(rep(0, 2 * np), nrow = 2 * np)
      arsumc <- matrix(rep(0, 2 * np * 2 * np), nrow = 2 * np)
      gi <- matrix(rep(0, 2 * np), nrow = 2 * np)
      arsumgfirstdev <- matrix(rep(0, 2 * np * np), nrow = 2 * np)
      firstdev <- matrix(rep(0, 2 * np * np), nrow = 2 * np)
      loc1 <- 0
      loc2 <- 0
      for (i in 1:nsub) {
        loc1 <- loc2 + 1
        loc2 <- loc1 + nobs[i] - 1
        yi <- as.matrix(Y[loc1:loc2, ])
        xi <- X[loc1:loc2, ]
        ni <- nrow(yi)
        cfixi<-cfix[loc1:loc2]
        m0 <- diag(ni)
        if(working=="ar1"){
          m1 <- matrix(rep(0, ni * ni), ni)
          for (k in 1:ni) {
            for (l in 1:ni) {
              if (abs(k - l) == 1)
                m1[k, l] <- 1
            }
          }
        }
        if(working=="exch"){
          m1 <- matrix(rep(1, ni * ni), ni) - m0
        }
        if(length(beta)!=1){
          ui <- mm(xi %*% beta+cfixi)
        }else{ ui <- mm(t(beta %*% t(xi))+cfixi)}
        fui <- log(ui)-log(1-ui)
        if(dim(ui)[[1]]==1)
        {
          vui <- diag(sqrt(1/ui)) %*% diag(sqrt(1/(1 - ui)))
          wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
          zi <- xi %*% fui_dev %*% vui %*% m1 %*% vui
        }
        if(dim(ui)[[1]]!=1){
          fui_dev <- diag((as.vector(ui))) %*% diag(as.vector(1 - ui))
          vui <- diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1 - ui))))
          wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
          zi <- t(xi) %*% fui_dev %*% vui %*% m1 %*% vui
        }
        gi0 <- (1/nsub) * wi %*% (yi - ui)
        gi1 <- (1/nsub) * zi %*% (yi - ui)
        gi[1:np, ] <- gi0
        gi[(np + 1):(2 * np), ] <- gi1
        arsumc <- arsumc + gi %*% t(gi)
        arsumg <- arsumg + gi
        di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
        di1 <- -(1/nsub) * zi %*% fui_dev %*% xi
        firstdev[1:np, ] <- di0
        firstdev[(np + 1):(2 * np), ] <- di1
        arsumgfirstdev <- arsumgfirstdev + firstdev
      }

      arcinv = geninv(arsumc)

      Q <- t(arsumg) %*% arcinv %*% arsumg #Qif
      #print(Q-QIF)
      QIF<-Q
      arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
      arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev
      pe<-as.numeric(pp_Lasso(beta,lambda_si,lambda_pl,T,dsi)/abs(beta))#scad penalty
      if(length(pe)!=1){
        sigma<- nsub*diag(pe)
      }else{sigma <- nsub*pe}
      #print("geninv(arqif2dev + sigma)")
      #print(arqif2dev + sigma)
      invarqif2dev <- geninv(arqif2dev + sigma)
      betanew <- beta - invarqif2dev %*% (arqif1dev + sigma %*% beta)
      a <- NULL
      b <- NULL
      for(i in 1:dim(betanew)[[1]]){ #cutoff
        #print("betanew")
        #print(betanew)
        if (abs(betanew[i])<1e-4) {
          betanew[i] <- 0 }
        else {
          a<-cbind(a,X[,i])
          b<-cbind(b,betanew[i])
        }
      }
      iid <- (1:length(BEST))[BEST!=0]
      BEST[iid] <- betanew
      betadiff <- t(betanew - beta)%*%(betanew - beta)
      # print("betadiff is ")
      #  print(betadiff)
      #print("betanew is ")
      #print(betanew)
      # print("QIF")
      #print(Q)
      if(sum(betanew)!=0){
        X <- as.matrix(a,nrow=nsub*ni)}
      np <- dim(X)[[2]]
      betanew <- as.numeric(b)
      iteration <- iteration + 1
      # print(iteration) shut off for now
    }
  }

  return(list(betanew=BEST,QIF=QIF))
}


# sqif_noApprox<-function(best,X,Y,splineest,lambda_si,lambda_pl,working,dsi,nsub,ni,T){
#   #best=b$coefficients[-1];X=XPQIF;Y=Y;lambda=lambda;working=corrQIF[corrStructure];nsub=n_sub;ni=n_obs
#   #best=best;X=XPQIF;Y=Y;lambda=0.2;working=corrQIF[corrStructure];nsub=n_sub;ni=n_obs
#   #best: initial estimate
#   #lambda : tuning parameter 'lambda'
#   #T : tuning parameter 'a'
#   #nsub : samle size
#   #ni : cluster size
#   #X : explanatory variable
#   #Y : response variable
#   #working : working correlation structure
#   X<-as.matrix(X)
#   Y<-as.matrix(Y)
#   np<-dim(X)[[2]]
#   id<-matrix(rep((1:nsub),each=ni),nsub*ni)
#   obs <- lapply(split(id, id), "length")
#   nobs <- as.numeric(obs) #number of cluster for each subject
#   nsub <- length(nobs) #number of subject
#   beta<-best#contained the intercept
#   betadiff <- 1
#   iteration <- 0
#   betanew <- beta
#   QIF=NULL
#   tol<-1e-10#very important this value have to be larger than simulation tol
#   maxiter<-1000
#   BEST <- as.matrix(rep(1,np),ncol=1)
#   QIF=0
#   if(working=="indep"){
#     while (betadiff > tol && iteration < maxiter && sum(betanew)!=0) {
#       beta <- betanew
#       sumg <- matrix(rep(0, np), nrow = np)
#       sumc <- matrix(rep(0, np * np), nrow = np)
#       arsumg <- matrix(rep(0, np), nrow = np)
#       arsumc <- matrix(rep(0, np * np), nrow = np)
#       gi <- matrix(rep(0, np), nrow = np)
#       arsumgfirstdev <- matrix(rep(0, np * np), nrow = np)
#       firstdev <- matrix(rep(0, np * np), nrow = np)
#       loc1 <- 0
#       loc2 <- 0
#
#       U=X%*%alphaold
#       bs=NULL
#       n=length(Y)
#       bsold = splineDesign(knots=knots, ord=(degree+1), x=U,derivs=rep(0, n), outer.ok=F)
#
#       for (i in 1:nsub) {
#         loc1 <- loc2 + 1
#         loc2 <- loc1 + nobs[i] - 1
#         yi <- as.matrix(Y[loc1:loc2, ])
#         xi <- X[loc1:loc2, ]
#         #cfixi<-cfix[loc1:loc2]
#         ni <- nrow(yi)
#         m0 <- diag(ni)
#         # print("cfix")
#         # print(cfix)
#
#         if(length(beta)!=1){
#           ui <-mm(bsold%*%splineest+zi%*%beta)
#          # ui <- mm(xi %*% beta)}#xi %*% beta+cfix
#         else{ ui <- mm(t(beta %*% t(xi)))}
#         fui <- log(ui)-log(1-ui) #link
#         if(dim(ui)[[1]]==1)
#         {
#           fui_dev <- diag(ui) %*% diag(1 - ui)
#           vui <- diag(sqrt(1/ui)) %*% diag(sqrt(1/(1 - ui)))
#           wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
#         } #fui_dev %*% vui %*% m0 %*% vui is identity maatrix
#         if(dim(ui)[[1]]!=1){
#           fui_dev <- diag((as.vector(ui))) %*% diag(as.vector(1 - ui))
#           vui <- diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1 - ui))))
#           wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
#         }
#         gi <- (1/nsub) * wi %*% (yi - ui)
#         arsumc <- arsumc + gi %*% t(gi)
#         arsumg <- arsumg + gi
#         di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
#         firstdev <- di0
#         arsumgfirstdev <- arsumgfirstdev + firstdev
#       } #end of looping through all subjects
#       #print(sprintf("det of weight matrix W %g",det(arsumc)))
#       print("Number of NAs in W inverse is ...")
#       print(sum(is.na(geninv(arsumc))))
#       arcinv = geninv(arsumc)
#       Q <- t(arsumg) %*% arcinv %*% arsumg
#       QIF<-Q
#       arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
#       arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev
#       pe<-as.numeric(pp_bin(beta,lambda_si,lambda_pl,T,dsi)/abs(beta))#scad penalty z,lambda_si,lambda_pl,a,dsi
#       if(length(pe)!=1){
#         sigma<- nsub*diag(pe)}
#       else{sigma <- nsub*pe}
#       #print("geninv(arqif2dev + sigma)")
#       #print(arqif2dev + sigma)
#       invarqif2dev <- geninv(arqif2dev + sigma)
#       betanew <- beta - invarqif2dev %*% (arqif1dev + sigma %*% beta)
#       a <- NULL
#       b <- NULL
#       for(i in 1:dim(betanew)[[1]]){
#         if (abs(betanew[i])<1e-4) {
#           betanew[i] <- 0 }
#         else {
#           a<-cbind(a,X[,i])
#           b<-cbind(b,betanew[i])
#         }
#       }
#       iid <- (1:length(BEST))[BEST!=0]
#       BEST[iid] <- betanew
#       betadiff <- t(betanew - beta)%*%(betanew - beta)
#       print(sprintf("betadiff is %g tol is %g",betadiff,tol))
#       print(betadiff)
#       print("betanew is ")
#       print(betanew)
#       if(sum(betanew)!=0){
#         X <- as.matrix(a,nrow=nsub*ni)}
#       np <- dim(X)[[2]]
#       betanew <- as.numeric(b)
#       iteration <- iteration + 1
#     }
#   }
#   if(working!="indep"){
#     while (betadiff > tol && iteration < maxiter && sum(betanew)!=0) {
#       beta <- betanew
#       sumg <- matrix(rep(0, 2 * np), nrow = 2 * np)
#       sumc <- matrix(rep(0, 2 * np * 2 * np), nrow = 2 * np)
#       arsumg <- matrix(rep(0, 2 * np), nrow = 2 * np)
#       arsumc <- matrix(rep(0, 2 * np * 2 * np), nrow = 2 * np)
#       gi <- matrix(rep(0, 2 * np), nrow = 2 * np)
#       arsumgfirstdev <- matrix(rep(0, 2 * np * np), nrow = 2 * np)
#       firstdev <- matrix(rep(0, 2 * np * np), nrow = 2 * np)
#       loc1 <- 0
#       loc2 <- 0
#       for (i in 1:nsub) {
#         loc1 <- loc2 + 1
#         loc2 <- loc1 + nobs[i] - 1
#         yi <- as.matrix(Y[loc1:loc2, ])
#         xi <- X[loc1:loc2, ]
#         ni <- nrow(yi)
#         #cfixi<-cfix[loc1:loc2]
#         m0 <- diag(ni)
#         if(working=="ar1"){
#           m1 <- matrix(rep(0, ni * ni), ni)
#           for (k in 1:ni) {
#             for (l in 1:ni) {
#               if (abs(k - l) == 1)
#                 m1[k, l] <- 1
#             }
#           }
#         }
#         if(working=="exch"){
#           m1 <- matrix(rep(1, ni * ni), ni) - m0
#         }
#         if(length(beta)!=1){
#           ui <-mm(bsold%*%splineest+zi%*%beta)
#         }else{ ui <- mm(t(beta %*% t(xi)))}
#         fui <- log(ui)-log(1-ui)
#         if(dim(ui)[[1]]==1)
#         {
#           vui <- diag(sqrt(1/ui)) %*% diag(sqrt(1/(1 - ui)))
#           wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
#           zi <- xi %*% fui_dev %*% vui %*% m1 %*% vui
#         }
#         if(dim(ui)[[1]]!=1){
#           fui_dev <- diag((as.vector(ui))) %*% diag(as.vector(1 - ui))
#           vui <- diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1 - ui))))
#           wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
#           zi <- t(xi) %*% fui_dev %*% vui %*% m1 %*% vui
#         }
#         gi0 <- (1/nsub) * wi %*% (yi - ui)
#         gi1 <- (1/nsub) * zi %*% (yi - ui)
#         # print("wi")
#         # print(wi)
#         # print("zi")
#         # print(zi)
#         gi[1:np, ] <- gi0
#         gi[(np + 1):(2 * np), ] <- gi1
#         arsumc <- arsumc + gi %*% t(gi)
#         arsumg <- arsumg + gi
#         di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
#         di1 <- -(1/nsub) * zi %*% fui_dev %*% xi
#         firstdev[1:np, ] <- di0
#         firstdev[(np + 1):(2 * np), ] <- di1
#         arsumgfirstdev <- arsumgfirstdev + firstdev
#       }
#       #print("***************gscore")
#       #print(dim(arsumg))
#       #print("dimbeta")
#       # print(length(betanew))
#       #  print(dim(arsumc))
#       #  arsumc=arsumc+0.0001*diag(nrow(arsumc))
#       # print(sprintf("det of weight matrix W %g",det(arsumc)))
#       # print("Number of NAs in W inverse is ...")
#       #  print(sum(is.na(geninv(arsumc))))
#       arcinv = geninv(arsumc)
#
#       Q <- t(arsumg) %*% arcinv %*% arsumg #Qif
#       #print(Q-QIF)
#       QIF<-Q
#       arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
#       arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev
#       pe<-as.numeric(pp_bin(beta,lambda_si,lambda_pl,T,dsi)/abs(beta))#scad penalty
#       if(length(pe)!=1){
#         sigma<- nsub*diag(pe)
#       }else{sigma <- nsub*pe}
#       #print("geninv(arqif2dev + sigma)")
#       #print(arqif2dev + sigma)
#       invarqif2dev <- geninv(arqif2dev + sigma)
#       betanew <- beta - invarqif2dev %*% (arqif1dev + sigma %*% beta)
#       a <- NULL
#       b <- NULL
#       for(i in 1:dim(betanew)[[1]]){ #cutoff
#         #print("betanew")
#         #print(betanew)
#         if (abs(betanew[i])<1e-4) {
#           betanew[i] <- 0 }
#         else {
#           a<-cbind(a,X[,i])
#           b<-cbind(b,betanew[i])
#         }
#       }
#       iid <- (1:length(BEST))[BEST!=0]
#       BEST[iid] <- betanew
#       betadiff <- t(betanew - beta)%*%(betanew - beta)
#       # print("betadiff is ")
#       #  print(betadiff)
#       #print("betanew is ")
#       #print(betanew)
#       # print("QIF")
#       #print(Q)
#       if(sum(betanew)!=0){
#         X <- as.matrix(a,nrow=nsub*ni)}
#       np <- dim(X)[[2]]
#       betanew <- as.numeric(b)
#       iteration <- iteration + 1
#       # print(iteration) shut off for now
#     }
#   }
#
#   return(list(betanew=BEST,QIF=QIF))
# }

##
######################################################
#Estimate parameter through our proposed method using#
#a linear approximation #
######################################################
scadplus <- function(beta,Y,X,n,J,p,working){
  y<-as.matrix(Y)
  x<-as.matrix(X)
  beta1<-matrix(beta, ncol=1)
  if(working=="indep"){
    gi <- matrix(rep(0, p), nrow = p)
    arsumg <- matrix(rep(0, p), nrow = p)
    arsumc <- matrix(rep(0, p * p), nrow = p)
    Yi<- matrix(rep(0, p), nrow = p)
    Xi<- matrix(rep(0, p * p), nrow = p)
    ys<- matrix(rep(0, p * n), ncol=1)
    xs<- matrix(rep(0, n * p * p), ncol=p)
    xss <-matrix(rep(0, n * p * p), ncol=p)
    m0 <- diag(J)
    for (i in 1:n)
    {
      yi<-y[(J*(i-1)+1):(J*i)]
      xi<-matrix(rep(1,(p*J)),J,p)
      xi<-x[(J*(i-1)+1):(J*i),]
      ui <- mm(xi %*% beta1)
      fui <- log(ui) - log(1 - ui)
      if(dim(ui)[[1]]==1)
      {
        fui_dev <- diag(ui) %*% diag(1 - ui)
        vui <- diag(sqrt(1/ui)) %*% diag(sqrt(1/(1 - ui)))
        wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
      }
      if(dim(ui)[[1]]!=1){
        fui_dev <- diag((as.vector(ui))) %*% diag(as.vector(1 - ui))
        vui <- diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1 - ui))))
        wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
      }
      gi <- (1/n) * wi %*% (yi - ui)
      ys[(p*(i-1)+1):(p*i)] <-(1/n)*(wi %*% t(t(xi) %*% fui_dev)%*% beta1)
      xs[(p*(i-1)+1):(p*i),] <-(1/n)* (wi %*% t(t(xi) %*% fui_dev))
      arsumc <- arsumc + gi %*% t(gi)
      arsumg <- arsumg + gi
    }
    arcinv = geninv(arsumc)
    #the square root of inverse of variance matrix
    A <- arcinv
    e <- eigen(A)
    V <- e$vectors
    #recovers A (up to rounding errors), and
    sqrtarcinv <- V %*% diag(sqrt(e$values)) %*% t(V)
    #find ystar and ystar
    for(i in 1:n)
    {
      ys[(p*(i-1)+1):(p*i)]<-sqrtarcinv%*%ys[(p*(i-1)+1):(p*i)]
      xs[(p*(i-1)+1):(p*i),]<-sqrtarcinv%*%xs[(p*(i-1)+1):(p*i),]
    }
    #final ystar (YS) and xdoublestar (xss)
    for (i in 1:p){
      xss[,i] <- xs[,i]
    }
    scad<-cbind(ys,xss) #new adjusted x and y by QIF method
    ys<-scad[,1]
    #SCAD
    second<-plus(xss,ys,m=3,intercept = FALSE)
    second.lasso<-plus(xss,ys,m=1,intercept = FALSE)
  }
  if(working!="indep"){
    gi <- matrix(rep(0, 2 * p), nrow = 2 * p)
    arsumg <- matrix(rep(0, 2 * p), nrow = 2 * p)
    arsumc <- matrix(rep(0, 2 * p * 2 * p), nrow = 2 * p)
    Yi<- matrix(rep(0, 2 * p), nrow = 2 * p)
    Xi<- matrix(rep(0, 2* p * p), nrow = 2 * p)
    ys<- matrix(rep(0, 2 * p * n), ncol=1)
    xs<- matrix(rep(0, n * 2 * p * p), ncol=p)
    xss <-matrix(rep(0, n * 2 * p * p), ncol=p)
    m0 <- diag(J)
    if(working=="ar1"){
      m1 <- matrix(rep(0, J * J), J)
      for (k in 1:J) {
        for (l in 1:J) {
          if (abs(k - l) == 1)
            m1[k, l] <- 1
        }
      }
    }
    if(working=="exch"){
      m1 <- matrix(rep(1, J * J), J) - m0
    }
    for (i in 1:n)
    {
      yi<-y[(J*(i-1)+1):(J*i)]
      xi<-matrix(rep(1,(p*J)),J,p)
      xi<-x[(J*(i-1)+1):(J*i),]
      ui <- mm(xi %*% beta1)
      fui <- log(ui) - log(1 - ui)
      if(dim(ui)[[1]]==1)
      {
        fui_dev <- diag(ui) %*% diag(1 - ui)
        vui <- diag(sqrt(1/ui)) %*% diag(sqrt(1/(1 - ui)))
        wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
        zi <- xi %*% fui_dev %*% vui %*% m1 %*% vui
      }
      if(dim(ui)[[1]]!=1){
        fui_dev <- diag((as.vector(ui))) %*% diag(as.vector(1 - ui))
        vui <- diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1 - ui))))
        wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
        zi <- t(xi) %*% fui_dev %*% vui %*% m1 %*% vui
      }
      gi0 <- (1/n) * wi %*% (yi - ui)
      gi1 <- (1/n) * zi %*% (yi - ui)
      gi[1:p, ] <- gi0
      gi[(p+1):(2*p), ] <- gi1
      Yi0 <-(1/n)*(wi %*% t(t(xi) %*% fui_dev)%*% beta1)
      Yi1 <-(1/n)*(zi %*% t(t(xi) %*% fui_dev)%*% beta1)
      ys[((2*p*i)-(2*p-1)):(2*p*i-p)] <- Yi0
      ys[((2*p*i)-(p-1)):(2*p*i)] <- Yi1
      Xi0 <-(1/n)* (wi %*% t(t(xi) %*% fui_dev))
      Xi1 <-(1/n)* (zi %*% t(t(xi) %*% fui_dev))
      xs[((2*p*i)-(2*p-1)):(2*p*i-p),] <- Xi0
      xs[((2*p*i)-(p-1)):(2*p*i),] <- Xi1
      arsumc <- arsumc + gi %*% t(gi)
      arsumg <- arsumg + gi
    }
    arcinv = geninv(arsumc)
    #the square root of inverse of variance matrix
    A <- arcinv
    e <- eigen(A)
    V <- e$vectors
    # recovers A (up to rounding errors), and
    sqrtarcinv <- V %*% diag(sqrt(e$values)) %*% t(V)
    #find ystar and ystar
    for(i in 1:n)
    {ys[((2*p*i)-(2*p-1)):(2*p*i)]<-sqrtarcinv%*%ys[((2*p*i)-(2*p-1)):(2*p*i)]
    xs[((2*p*i)-(2*p-1)):(2*p*i),]<-sqrtarcinv%*%xs[((2*p*i)-(2*p-1)):(2*p*i),]
    }
    #final ystar (YS) and xdoublestar (xss)
    for (i in 1:p){
      xss[,i] <- xs[,i]
    }
    scad<-cbind(ys,xss) #new adjusted x and y by QIF method
    ys<-scad[,1]
    #SCAD
    second<-plus(xss,ys,m=3,intercept = FALSE)
    second.lasso<-plus(xss,ys,m=1,intercept = FALSE)
  }
  #Find the optimal tuning parameter lambda
  extracted.coef <- coef(second, lam = sort(second$lam.path,decreasing=TRUE))
  df<-second$dim
  y<-as.matrix(Y)
  x<-as.matrix(X)
  k<-nrow(extracted.coef)
  BIQIF <- c(rep(0,k))
  Q <- c(rep(0,k))
  if(working=="indep"){
    for(j in 1:k)
    {
      beta1<-matrix(extracted.coef[j,], nrow=p)
      gi <- matrix(rep(0, p), nrow = p)
      arsumg <- matrix(rep(0, p), nrow = p)
      arsumc <- matrix(rep(0, p * p), nrow = p)
      for (i in 1:n)
      {
        yi<-y[(J*(i-1)+1):(J*i)]
        xi<-matrix(rep(1,(p*J)),J,p)
        xi<-x[(J*(i-1)+1):(J*i),]
        ui <- mm(xi %*% beta1)
        fui <- log(ui) - log(1 - ui)
        fui_dev <- diag(as.vector(ui)) %*% diag(as.vector(1 - ui))
        vui <- diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1-ui))))
        wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
        gi <- (1/n) * wi %*% (yi - ui)
        arsumc <- arsumc + gi %*% t(gi)
        arsumg <- arsumg + gi
      }
      arcinv = geninv(arsumc)
      Q[j] <- t(arsumg) %*% arcinv %*% arsumg
      BIQIF[j]<-Q[j]+df[j]*log(n)
    }
  }
  if(working!="indep"){
    for(j in 1:k)
    {
      beta1<-matrix(extracted.coef[j,], nrow=p)
      gi <- matrix(rep(0, 2 * p), nrow = 2*p)
      arsumg <- matrix(rep(0, 2*p), nrow = 2*p)
      arsumc <- matrix(rep(0, 2*p * 2*p), nrow = 2*p)
      for (i in 1:n)
      {
        yi<-y[(J*(i-1)+1):(J*i)]
        xi<-matrix(rep(1,(p*J)),J,p)
        xi<-x[(J*(i-1)+1):(J*i),]
        ui <- mm(xi %*% beta1)
        fui <- log(ui) - log(1 - ui)
        fui_dev <- diag(as.vector(ui)) %*% diag(as.vector(1 - ui))
        vui <- diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1-ui))))
        wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
        zi <- t(xi) %*% fui_dev %*% vui %*% m1 %*% vui
        gi0 <- (1/n) * wi %*% (yi - ui)
        gi1 <- (1/n) * zi %*% (yi - ui)
        gi[1:p, ] <- gi0
        gi[(p+1):(2*p), ] <- gi1
        arsumc <- arsumc + gi %*% t(gi)
        arsumg <- arsumg + gi
      }
      arcinv <- geninv(arsumc)
      Q[j] <- t(arsumg) %*% arcinv %*% arsumg
      BIQIF[j]<-Q[j]+df[j]*log(n)
    }
  }
  #Find the best estimate by minimizing BIQIF
  mmm <- (1:k)[BIQIF==min(BIQIF)]
  BEST <- extracted.coef[mmm[1],]
  return(list(BEST=BEST))
}
######################################################
#Estimate parameter through our proposed method using#
#a local quadratic approximation for continuous response #
######################################################

gausqifUnbalanced<-function(best,X,Y,id,lambda,working,nsub,pindex,T){
  #best=initials;X=XLin;lambda=Lingrid[l];working=corrQIF[corrStructure];nsub=n_sub;pindex=1
  #best=initials;X=XZ;Y=Y;lambda=bestLinlambda[i];working=corrQIF[corrStructure];nsub=n_sub;ni=n_obs
  #best: initial estimate
  #lambda : tuning parameter 'lambda'
  #T : tuning parameter 'a'
  #nsub : samle size
  #ni : cluster size
  #X : explanatory variable
  #Y : response variable
  #working : working correlation structure best=c(true.alpha,true.beta)
  X<-as.matrix(X)
  Y<-as.matrix(Y)
  np<-dim(X)[[2]]
  #id<-matrix(rep((1:nsub),each=ni),nsub*ni)
  obs <- lapply(split(id, id), "length")
  nobs <- as.numeric(obs) #number of cluster for each subject
  nsub <- length(nobs) #number of subject
  beta<-best#contained the intercept
  betadiff <- 1
  iteration <- 0
  betanew <- beta
  tol<-1e-10#very important this value have to be larger than simulation tol
  maxiter<-1000
  BEST <- as.matrix(rep(1,np),ncol=1)
  if(working=="indep"){
    while (betadiff > tol && iteration < maxiter && sum(betanew)!=0) {
      beta <- betanew
      sumg <- matrix(rep(0, np), nrow = np)
      sumc <- matrix(rep(0, np * np), nrow = np)
      arsumg <- matrix(rep(0, np), nrow = np)
      arsumc <- matrix(rep(0, np * np), nrow = np)
      gi <- matrix(rep(0, np), nrow = np)
      arsumgfirstdev <- matrix(rep(0, np * np), nrow = np)
      firstdev <- matrix(rep(0, np * np), nrow = np)
      loc1 <- 0
      loc2 <- 0
      for (i in 1:nsub) {
        ni=nobs[i]
        loc1 <- loc2 + 1
        loc2 <- loc1 + nobs[i] - 1
        yi <- as.matrix(Y[loc1:loc2, ])
        xi <- X[loc1:loc2, ]
        ni <- nrow(yi)
        m0 <- diag(ni)
        if(length(beta)!=1){
          ui <- xi %*% beta}
        else{ ui <- t(beta %*% t(xi))}
        fui <- ui
        if(dim(ui)[[1]]==1)
        {
          fui_dev <- diag(ni)#diag(ui) %*% diag(1 - ui)
          vui <- diag(ni)#diag(sqrt(1/ui)) %*% diag(sqrt(1/(1 - ui)))
          # wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
          wi<-xi %*% m0
        } #fui_dev %*% vui %*% m0 %*% vui is identity maatrix
        if(dim(ui)[[1]]!=1){
          fui_dev <- diag(ni)#diag((as.vector(ui))) %*% diag(as.vector(1 - ui))
          vui <- diag(ni)#diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1 - ui))))
          # wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
          wi <- t(xi) %*%  m0
        }
        gi <- (1/nsub) * wi %*% (yi - ui)
        arsumc <- arsumc + gi %*% t(gi)
        arsumg <- arsumg +gi
        #di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
        di0 <- -(1/nsub) * wi  %*% xi
        firstdev <- di0
        arsumgfirstdev <- arsumgfirstdev + firstdev
      }
      arcinv = geninv(arsumc)
      Q <- t(arsumg) %*% arcinv %*% arsumg
      arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
      arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev
      pe<-as.numeric(pp_1lambda(beta,lambda,T)/abs(beta))#scad penalty
      if(is.null(pindex)!=TRUE) {
        pe[pindex]<-0
      }
      if(length(pe)!=1){
        sigma<- nsub*diag(pe)}
      else{sigma <- nsub*pe}
      #print("geninv(arqif2dev + sigma)")
      #print(arqif2dev + sigma)
      invarqif2dev <- geninv(arqif2dev + sigma)
      betanew <- beta - invarqif2dev %*% (arqif1dev + sigma %*% beta)
      a <- NULL
      b <- NULL
      for(i in 1:dim(betanew)[[1]]){
        if (abs(betanew[i])<1e-4) {
          betanew[i] <- 0 }
        else {
          a<-cbind(a,X[,i])
          b<-cbind(b,betanew[i])
        }
      }
      iid <- (1:length(BEST))[BEST!=0]
      BEST[iid] <- betanew
      betadiff <- t(betanew - beta)%*%(betanew - beta)
      if(sum(betanew)!=0){
        X <- as.matrix(a,nrow=nsub*ni)}
      np <- dim(X)[[2]]
      betanew <- as.numeric(b)
      iteration <- iteration + 1
    }
  }
  if(working!="indep"){
    while (betadiff > tol && iteration < maxiter && sum(betanew)!=0) {
      beta <- betanew
      sumg <- matrix(rep(0, 2 * np), nrow = 2 * np)
      sumc <- matrix(rep(0, 2 * np * 2 * np), nrow = 2 * np)
      arsumg <- matrix(rep(0, 2 * np), nrow = 2 * np)
      arsumc <- matrix(rep(0, 2 * np * 2 * np), nrow = 2 * np)
      gi <- matrix(rep(0, 2 * np), nrow = 2 * np)
      arsumgfirstdev <- matrix(rep(0, 2 * np * np), nrow = 2 * np)
      firstdev <- matrix(rep(0, 2 * np * np), nrow = 2 * np)
      loc1 <- 0
      loc2 <- 0
      for (i in 1:nsub) {#i=1
        ni=nobs[i]
        loc1 <- loc2 + 1
        loc2 <- loc1 + nobs[i] - 1
        yi <- as.matrix(Y[loc1:loc2, ])
        xi <- X[loc1:loc2, ]
        ni <- nrow(yi)
        m0 <- diag(ni)
        if(working=="ar1"){
          m1 <- matrix(rep(0, ni * ni), ni)
          for (k in 1:ni) {
            for (l in 1:ni) {
              if (abs(k - l) == 1)
                m1[k, l] <- 1
            }
          }
        }
        if(working=="exch"){
          m1 <- matrix(rep(1, ni * ni), ni) - m0
        }
        #first compute residual
        if(length(beta)!=1){
          ui <- xi %*% beta #mm(xi %*% beta)
        }else{ ui <- t(beta %*% t(xi)) }#mm(t(beta %*% t(xi)))}

        fui <- ui#log(ui)-log(1-ui)
        if(dim(ui)[[1]]==1)
        {
          fui_dev <- diag(nrow(ui))#diag(ui) %*% diag(1 - ui)
          vui <- diag(nrow(ui))#diag(sqrt(1/ui)) %*% diag(sqrt(1/(1 - ui)))
          wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
          zi <- xi %*% fui_dev %*% vui %*% m1 %*% vui
          # wi <- xi %*% m0
          #  zi <- xi %*% m1


        }
        if(dim(ui)[[1]]!=1){

          fui_dev <-diag(nrow(ui)) #binary->diag((as.vector(ui))) %*% diag(as.vector(1 - ui))
          vui <- diag(nrow(ui)) #diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1 - ui)))) #A-1/2
          wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui #first row for gscore D'A^1/2 Mo A^1/2 with D'=X'A
          zi <- t(xi) %*% fui_dev %*% vui %*% m1 %*% vui #second row for gscore D'A^1/2 M1 A^1/2
          #wi <- t(xi)%*% m0
          #zi <- t(xi)%*% m1

        }
        #building the g score
        gi0 <- (1/nsub) * wi %*% (yi - ui)
        gi1 <- (1/nsub) * zi %*% (yi - ui)
        gi[1:np, ] <- gi0
        gi[(np + 1):(2 * np), ] <- gi1
        #sum W^-1 over all subjects
        arsumc <- arsumc + gi %*% t(gi)
        #sum gscores over all subjects
        arsumg <- arsumg + gi
        # di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
        #di1 <- -(1/nsub) * zi %*% fui_dev %*% xi
        di0 <- -(1/nsub) * wi  %*% xi
        di1 <- -(1/nsub) * zi %*% xi
        firstdev[1:np, ] <- di0
        firstdev[(np + 1):(2 * np), ] <- di1
        arsumgfirstdev <- arsumgfirstdev + firstdev
      }
      #print(arsumc)
      #print(dim(arsumc))
      #arsumc=arsumc+0.0001*diag(nrow(arcsum))
      #print("arsumc")
      arcinv = geninv(arsumc)
      #print(arcinv)
      #calculating QIF
      Q <- t(arsumg) %*% arcinv %*% arsumg #Qif
      arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
      arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev

      pe<-as.numeric(pp_1lambda(beta,lambda,T)/abs(beta))#scad penalty
      if(is.null(pindex)!=TRUE) {
        pe[pindex]<-0
      }

      if(length(pe)!=1){
        sigma<- nsub*diag(pe)
      }else{sigma <- nsub*pe}
      invarqif2dev <- geninv(arqif2dev + sigma)#,tol=0)
      betanew <- beta - invarqif2dev %*% (arqif1dev + sigma %*% beta)
      a <- NULL
      b <- NULL
      for(i in 1:dim(betanew)[[1]]){
        if (abs(betanew[i])<1e-4) {
          betanew[i] <- 0 }
        else {
          a<-cbind(a,X[,i])
          b<-cbind(b,betanew[i])
        }
      }
      iid <- (1:length(BEST))[BEST!=0]
      BEST[iid] <- betanew
      betadiff <- t(betanew - beta)%*%(betanew - beta)
      if(sum(betanew)!=0){
        X <- as.matrix(a,nrow=nsub*ni)}
      np <- dim(X)[[2]]
      betanew <- as.numeric(b)
      iteration <- iteration + 1
      # print(iteration) shut off for now
    }
  }
  return(list(betanew=BEST,QIF=Q))
}

gausqifLASSO<-function(best,X,Y,lambda,working,nsub,ni,T){
  #best=initials;X=XZ;Y=Y;lambda=bestLinlambda[i];working=corrQIF[corrStructure];nsub=n_sub;ni=n_obs
  #best: initial estimate
  #lambda : tuning parameter 'lambda'
  #T : tuning parameter 'a'
  #nsub : samle size
  #ni : cluster size
  #X : explanatory variable
  #Y : response variable
  #working : working correlation structure best=c(true.alpha,true.beta)
  X<-as.matrix(X)
  Y<-as.matrix(Y)
  np<-dim(X)[[2]]
  id<-matrix(rep((1:nsub),each=ni),nsub*ni)
  obs <- lapply(split(id, id), "length")
  nobs <- as.numeric(obs) #number of cluster for each subject
  nsub <- length(nobs) #number of subject
  beta<-best#contained the intercept
  betadiff <- 1
  iteration <- 0
  betanew <- beta
  tol<-1e-10#very important this value have to be larger than simulation tol
  maxiter<-1000
  BEST <- as.matrix(rep(1,np),ncol=1)
  if(working=="indep"){
    while (betadiff > tol && iteration < maxiter && sum(betanew)!=0) {
      beta <- betanew
      sumg <- matrix(rep(0, np), nrow = np)
      sumc <- matrix(rep(0, np * np), nrow = np)
      arsumg <- matrix(rep(0, np), nrow = np)
      arsumc <- matrix(rep(0, np * np), nrow = np)
      gi <- matrix(rep(0, np), nrow = np)
      arsumgfirstdev <- matrix(rep(0, np * np), nrow = np)
      firstdev <- matrix(rep(0, np * np), nrow = np)
      loc1 <- 0
      loc2 <- 0
      for (i in 1:nsub) {
        loc1 <- loc2 + 1
        loc2 <- loc1 + nobs[i] - 1
        yi <- as.matrix(Y[loc1:loc2, ])
        xi <- X[loc1:loc2, ]
        ni <- nrow(yi)
        m0 <- diag(ni)
        if(length(beta)!=1){
          ui <- xi %*% beta}
        else{ ui <- t(beta %*% t(xi))}
        fui <- ui
        if(dim(ui)[[1]]==1)
        {
          fui_dev <- diag(nrow(ui))#diag(ui) %*% diag(1 - ui)
          vui <- diag(nrow(ui))#diag(sqrt(1/ui)) %*% diag(sqrt(1/(1 - ui)))
          # wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
          wi<-xi %*% m0
        } #fui_dev %*% vui %*% m0 %*% vui is identity maatrix
        if(dim(ui)[[1]]!=1){
          fui_dev <- diag(nrow(ui))#diag((as.vector(ui))) %*% diag(as.vector(1 - ui))
          vui <- diag(nrow(ui))#diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1 - ui))))
          # wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
          wi <- t(xi) %*%  m0
        }
        gi <- (1/nsub) * wi %*% (yi - ui)
        arsumc <- arsumc + gi %*% t(gi)
        arsumg <- arsumg +gi
        #di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
        di0 <- -(1/nsub) * wi  %*% xi
        firstdev <- di0
        arsumgfirstdev <- arsumgfirstdev + firstdev
      }
      arcinv = geninv(arsumc)
      Q <- t(arsumg) %*% arcinv %*% arsumg
      arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
      arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev
      pe<-as.numeric(pplasso(beta,lambda,T)/abs(beta))#scad penalty
      if(length(pe)!=1){
        sigma<- nsub*diag(pe)}
      else{sigma <- nsub*pe}
      #print("geninv(arqif2dev + sigma)")
      #print(arqif2dev + sigma)
      invarqif2dev <- geninv(arqif2dev + sigma)
      betanew <- beta - invarqif2dev %*% (arqif1dev + sigma %*% beta)
      a <- NULL
      b <- NULL
      for(i in 1:dim(betanew)[[1]]){
        if (abs(betanew[i])<1e-4) {
          betanew[i] <- 0 }
        else {
          a<-cbind(a,X[,i])
          b<-cbind(b,betanew[i])
        }
      }
      iid <- (1:length(BEST))[BEST!=0]
      BEST[iid] <- betanew
      betadiff <- t(betanew - beta)%*%(betanew - beta)
      if(sum(betanew)!=0){
        X <- as.matrix(a,nrow=nsub*ni)}
      np <- dim(X)[[2]]
      betanew <- as.numeric(b)
      iteration <- iteration + 1
    }
  }
  if(working!="indep"){
    while (betadiff > tol && iteration < maxiter && sum(betanew)!=0) {
      beta <- betanew
      sumg <- matrix(rep(0, 2 * np), nrow = 2 * np)
      sumc <- matrix(rep(0, 2 * np * 2 * np), nrow = 2 * np)
      arsumg <- matrix(rep(0, 2 * np), nrow = 2 * np)
      arsumc <- matrix(rep(0, 2 * np * 2 * np), nrow = 2 * np)
      gi <- matrix(rep(0, 2 * np), nrow = 2 * np)
      arsumgfirstdev <- matrix(rep(0, 2 * np * np), nrow = 2 * np)
      firstdev <- matrix(rep(0, 2 * np * np), nrow = 2 * np)
      loc1 <- 0
      loc2 <- 0
      for (i in 1:nsub) {#i=1
        loc1 <- loc2 + 1
        loc2 <- loc1 + nobs[i] - 1
        yi <- as.matrix(Y[loc1:loc2, ])
        xi <- X[loc1:loc2, ]
        ni <- nrow(yi)
        m0 <- diag(ni)
        if(working=="ar1"){
          m1 <- matrix(rep(0, ni * ni), ni)
          for (k in 1:ni) {
            for (l in 1:ni) {
              if (abs(k - l) == 1)
                m1[k, l] <- 1
            }
          }
        }
        if(working=="exch"){
          m1 <- matrix(rep(1, ni * ni), ni) - m0
        }
        #first compute residual
        if(length(beta)!=1){
          ui <- xi %*% beta #mm(xi %*% beta)
        }else{ ui <- t(beta %*% t(xi)) }#mm(t(beta %*% t(xi)))}

        fui <- ui#log(ui)-log(1-ui)
        if(dim(ui)[[1]]==1)
        {
          fui_dev <- diag(nrow(ui))#diag(ui) %*% diag(1 - ui)
          vui <- diag(nrow(ui))#diag(sqrt(1/ui)) %*% diag(sqrt(1/(1 - ui)))
          wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
          zi <- xi %*% fui_dev %*% vui %*% m1 %*% vui
          # wi <- xi %*% m0
          #  zi <- xi %*% m1


        }
        if(dim(ui)[[1]]!=1){

          fui_dev <-diag(nrow(ui)) #binary->diag((as.vector(ui))) %*% diag(as.vector(1 - ui))
          vui <- diag(nrow(ui)) #diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1 - ui)))) #A-1/2
          wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui #first row for gscore D'A^1/2 Mo A^1/2 with D'=X'A
          zi <- t(xi) %*% fui_dev %*% vui %*% m1 %*% vui #second row for gscore D'A^1/2 M1 A^1/2
          #wi <- t(xi)%*% m0
          #zi <- t(xi)%*% m1

        }
        #building the g score
        gi0 <- (1/nsub) * wi %*% (yi - ui)
        gi1 <- (1/nsub) * zi %*% (yi - ui)
        gi[1:np, ] <- gi0
        gi[(np + 1):(2 * np), ] <- gi1
        #sum W^-1 over all subjects
        arsumc <- arsumc + gi %*% t(gi)
        #sum gscores over all subjects
        arsumg <- arsumg + gi
        # di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
        #di1 <- -(1/nsub) * zi %*% fui_dev %*% xi
        di0 <- -(1/nsub) * wi  %*% xi
        di1 <- -(1/nsub) * zi %*% xi
        firstdev[1:np, ] <- di0
        firstdev[(np + 1):(2 * np), ] <- di1
        arsumgfirstdev <- arsumgfirstdev + firstdev
      }
      #print(arsumc)
      #print(dim(arsumc))
      #arsumc=arsumc+0.0001*diag(nrow(arcsum))
      #print("arsumc")
      arcinv = geninv(arsumc)
      #print(arcinv)
      #calculating QIF
      Q <- t(arsumg) %*% arcinv %*% arsumg #Qif
      arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
      arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev
      pe<-as.numeric(pplasso(beta,lambda,T)/abs(beta))#scad penalty
      if(length(pe)!=1){
        sigma<- nsub*diag(pe)
      }else{sigma <- nsub*pe}
      # print("arqif2dev + sigma for geninv(arqif2dev + sigma)")
      #print(arqif2dev + sigma)
      invarqif2dev <- geninv(arqif2dev + sigma)#,tol=0)
      betanew <- beta - invarqif2dev %*% (arqif1dev + sigma %*% beta)
      a <- NULL
      b <- NULL
      for(i in 1:dim(betanew)[[1]]){
        if (abs(betanew[i])<1e-4) {
          betanew[i] <- 0 }
        else {
          a<-cbind(a,X[,i])
          b<-cbind(b,betanew[i])
        }
      }
      iid <- (1:length(BEST))[BEST!=0]
      BEST[iid] <- betanew
      betadiff <- t(betanew - beta)%*%(betanew - beta)
      if(sum(betanew)!=0){
        X <- as.matrix(a,nrow=nsub*ni)}
      np <- dim(X)[[2]]
      betanew <- as.numeric(b)
      iteration <- iteration + 1
      # print(iteration) shut off for now
    }
  }
  return(list(betanew=BEST,QIF=Q))
}


gausqif<-function(best,X,Y,lambda_si,lambda_pl,working,dsi,nsub,ni,T){
  #best=initials;X=XZ;Y=Y;lambda=bestLinlambda[i];working=corrQIF[corrStructure];nsub=n_sub;ni=n_obs
  #best: initial estimate
  #lambda : tuning parameter 'lambda'
  #T : tuning parameter 'a'
  #nsub : samle size
  #ni : cluster size
  #X : explanatory variable
  #Y : response variable
  #working : working correlation structure best=c(true.alpha,true.beta)
  X<-as.matrix(X)
  Y<-as.matrix(Y)
  np<-dim(X)[[2]]
  id<-matrix(rep((1:nsub),each=ni),nsub*ni)
  obs <- lapply(split(id, id), "length")
  nobs <- as.numeric(obs) #number of cluster for each subject
  nsub <- length(nobs) #number of subject
  beta<-best#contained the intercept
  betadiff <- 1
  iteration <- 0
  betanew <- beta
  tol<-1e-10#very important this value have to be larger than simulation tol
  maxiter<-1000
  BEST <- as.matrix(rep(1,np),ncol=1)
  if(working=="indep"){
    while (betadiff > tol && iteration < maxiter && sum(betanew)!=0) {
      beta <- betanew
      sumg <- matrix(rep(0, np), nrow = np)
      sumc <- matrix(rep(0, np * np), nrow = np)
      arsumg <- matrix(rep(0, np), nrow = np)
      arsumc <- matrix(rep(0, np * np), nrow = np)
      gi <- matrix(rep(0, np), nrow = np)
      arsumgfirstdev <- matrix(rep(0, np * np), nrow = np)
      firstdev <- matrix(rep(0, np * np), nrow = np)
      loc1 <- 0
      loc2 <- 0
      for (i in 1:nsub) {
        loc1 <- loc2 + 1
        loc2 <- loc1 + nobs[i] - 1
        yi <- as.matrix(Y[loc1:loc2, ])
        xi <- X[loc1:loc2, ]
        ni <- nrow(yi)
        m0 <- diag(ni)
        if(length(beta)!=1){
          ui <- xi %*% beta}
        else{ ui <- t(beta %*% t(xi))}
        fui <- ui
        if(dim(ui)[[1]]==1)
        {
          fui_dev <- diag(ni)#diag(ui) %*% diag(1 - ui)
          vui <- diag(ni)#diag(sqrt(1/ui)) %*% diag(sqrt(1/(1 - ui)))
          # wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
          wi<-xi %*% m0
        } #fui_dev %*% vui %*% m0 %*% vui is identity maatrix
        if(dim(ui)[[1]]!=1){
          fui_dev <- diag(ni)#diag((as.vector(ui))) %*% diag(as.vector(1 - ui))
          vui <- diag(ni)#diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1 - ui))))
          # wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
          wi <- t(xi) %*%  m0
        }
        gi <- (1/nsub) * wi %*% (yi - ui)
        arsumc <- arsumc + gi %*% t(gi)
        arsumg <- arsumg +gi
        #di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
        di0 <- -(1/nsub) * wi  %*% xi
        firstdev <- di0
        arsumgfirstdev <- arsumgfirstdev + firstdev
      }
      arcinv = geninv(arsumc)
      #print("befor inverse is ")
      #print(arsumc)
      # print("inverse is ")
      #print(arcinv)
      Q <- t(arsumg) %*% arcinv %*% arsumg
      arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
      arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev
      pe<-as.numeric(pp(beta,lambda_si,lambda_pl,T,dsi)/abs(beta))#scad penalty
      #print(pe)
      if(length(pe)!=1){
        sigma<- nsub*diag(pe)}
      else{sigma <- nsub*pe}
      #print("geninv(arqif2dev + sigma)")
      #print(arqif2dev + sigma)
      invarqif2dev <- geninv(arqif2dev + sigma)
      betanew <- beta - invarqif2dev %*% (arqif1dev + sigma %*% beta)
      a <- NULL
      b <- NULL
      for(i in 1:dim(betanew)[[1]]){
        if (abs(betanew[i])<1e-4) {
          betanew[i] <- 0 }
        else {
          a<-cbind(a,X[,i])
          b<-cbind(b,betanew[i])
        }
      }
      iid <- (1:length(BEST))[BEST!=0]
      BEST[iid] <- betanew
      betadiff <- t(betanew - beta)%*%(betanew - beta)
      if(sum(betanew)!=0){
        X <- as.matrix(a,nrow=nsub*ni)}
      np <- dim(X)[[2]]
      betanew <- as.numeric(b)
      iteration <- iteration + 1
    }
  }
  if(working!="indep"){
    while (betadiff > tol && iteration < maxiter && sum(betanew)!=0) {
      beta <- betanew
      sumg <- matrix(rep(0, 2 * np), nrow = 2 * np)
      sumc <- matrix(rep(0, 2 * np * 2 * np), nrow = 2 * np)
      arsumg <- matrix(rep(0, 2 * np), nrow = 2 * np)
      arsumc <- matrix(rep(0, 2 * np * 2 * np), nrow = 2 * np)
      gi <- matrix(rep(0, 2 * np), nrow = 2 * np)
      arsumgfirstdev <- matrix(rep(0, 2 * np * np), nrow = 2 * np)
      firstdev <- matrix(rep(0, 2 * np * np), nrow = 2 * np)
      loc1 <- 0
      loc2 <- 0
      for (i in 1:nsub) {#i=1
        loc1 <- loc2 + 1
        loc2 <- loc1 + nobs[i] - 1
        yi <- as.matrix(Y[loc1:loc2, ])
        xi <- X[loc1:loc2, ]
        ni <- nrow(yi)
        m0 <- diag(ni)
        if(working=="ar1"){
          m1 <- matrix(rep(0, ni * ni), ni)
          for (k in 1:ni) {
            for (l in 1:ni) {
              if (abs(k - l) == 1)
                m1[k, l] <- 1
            }
          }
        }
        if(working=="exch"){
          m1 <- matrix(rep(1, ni * ni), ni) - m0
        }
        #first compute residual
        if(length(beta)!=1){
          ui <- xi %*% beta #mm(xi %*% beta)
        }else{ ui <- t(beta %*% t(xi)) }#mm(t(beta %*% t(xi)))}

        fui <- ui#log(ui)-log(1-ui)
        if(dim(ui)[[1]]==1)
        {
          fui_dev <- diag(nrow(ui))#diag(ui) %*% diag(1 - ui)
          vui <- diag(nrow(ui))#diag(sqrt(1/ui)) %*% diag(sqrt(1/(1 - ui)))
          wi <- xi %*% fui_dev %*% vui %*% m0 %*% vui
          zi <- xi %*% fui_dev %*% vui %*% m1 %*% vui
          # wi <- xi %*% m0
          #  zi <- xi %*% m1


        }
        if(dim(ui)[[1]]!=1){

          fui_dev <-diag(nrow(ui)) #binary->diag((as.vector(ui))) %*% diag(as.vector(1 - ui))
          vui <- diag(nrow(ui)) #diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1 - ui)))) #A-1/2
          wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui #first row for gscore D'A^1/2 Mo A^1/2 with D'=X'A
          zi <- t(xi) %*% fui_dev %*% vui %*% m1 %*% vui #second row for gscore D'A^1/2 M1 A^1/2
          #wi <- t(xi)%*% m0
          #zi <- t(xi)%*% m1

        }
        #building the g score
        gi0 <- (1/nsub) * wi %*% (yi - ui)
        gi1 <- (1/nsub) * zi %*% (yi - ui)
        gi[1:np, ] <- gi0
        gi[(np + 1):(2 * np), ] <- gi1
        #sum W^-1 over all subjects
        arsumc <- arsumc + gi %*% t(gi)
        #sum gscores over all subjects
        arsumg <- arsumg + gi
        # di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
        #di1 <- -(1/nsub) * zi %*% fui_dev %*% xi
        di0 <- -(1/nsub) * wi  %*% xi
        di1 <- -(1/nsub) * zi %*% xi
        firstdev[1:np, ] <- di0
        firstdev[(np + 1):(2 * np), ] <- di1
        arsumgfirstdev <- arsumgfirstdev + firstdev
      }
      #print(arsumc)
      #print(dim(arsumc))
      #arsumc=arsumc+0.0001*diag(nrow(arcsum))
      #print("arsumc")
      #  print("before inverse is ")
      #  print(arsumc)
      # print("inverse is ")
      #  print(arcinv)
      arcinv = geninv(arsumc)
      #print(arcinv)
      #calculating QIF
      Q <- t(arsumg) %*% arcinv %*% arsumg #Qif
      arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
      arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev
      pe<-as.numeric(pp(beta,lambda_si,lambda_pl,T,dsi)/abs(beta))#scad penalty
      # print(pe)
      if(length(pe)!=1){
        sigma<- nsub*diag(pe)
      }else{sigma <- nsub*pe}
      invarqif2dev <- geninv(arqif2dev + sigma)#,tol=0)
      betanew <- beta - invarqif2dev %*% (arqif1dev + sigma %*% beta)
      a <- NULL
      b <- NULL
      for(i in 1:dim(betanew)[[1]]){
        if (abs(betanew[i])<1e-4) {
          betanew[i] <- 0 }
        else {
          a<-cbind(a,X[,i])
          b<-cbind(b,betanew[i])
        }
      }
      iid <- (1:length(BEST))[BEST!=0]
      BEST[iid] <- betanew
      betadiff <- t(betanew - beta)%*%(betanew - beta)
      if(sum(betanew)!=0){
        X <- as.matrix(a,nrow=nsub*ni)}
      np <- dim(X)[[2]]
      betanew <- as.numeric(b)
      iteration <- iteration + 1
      # print(iteration) shut off for now
    }
  }
  return(list(betanew=BEST,QIF=Q))
}



