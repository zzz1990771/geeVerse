#main method
#Quantile Longitudinal GEE SCAD
#####x,y data
#####betaint: intial estimate for beta
#####nk: the number of observations in each subject
#####worktype: working matrix
#####f0: estimate for the density of error
qlingee_scad<-function(x,y,betaint,nk,worktype,f0=rep(1,length(y)),tau=0.5,lambda=0.1,max_it = 100,cutoff=10^-1){
  x_all = x
  cn = c(0, cumsum(nk));nsub=length(nk);nx=dim(x)[2];
  N=sum(nk);
  betaold=beta_all=beta=betaint;
  ul=x%*%beta;
  ghat=ul
  iteration = 0
  w = 1
  diff2a=1
  removed_ind = c()

  kk=0 #maybe remove

  while (max(abs(diff2a))>0.001 & iteration<max_it){
    iteration=iteration+1

    new_removed_ind <- which(abs(beta)<=cutoff)
    if(length(new_removed_ind)!=0 &&
       length(c(removed_ind,new_removed_ind))< -1+length(beta_all)){
      #update removed_ind
      if(is.null(removed_ind)){
        removed_ind = new_removed_ind
      }else{
        removed_ind = c(removed_ind,
                        (1:length(beta_all))[-removed_ind][new_removed_ind])
      }
      #update X and beta
      if(!is.null(removed_ind)){
        x = x_all[,-removed_ind]
        beta_all[removed_ind] = 0
        beta = beta_all[-removed_ind]
        beta_all[-removed_ind] = beta
      }else{
        x = x_all
        beta = beta_all
        beta_all = beta
      }


    }

    resid<-NULL;R<-sumuv<-matrix(0,max(nk),max(nk));

    resid=0.5-(y-ghat<0);######## used to calculate the correlation matrix

    sd_resid=resid/sqrt(tau-tau^2) ##########standard the resid

    #sd_resida=NULL;for (i in 1:nsub) sd_resida[[i]]=sd_resid[(cn[i]+1):cn[i+1]]


    if (worktype=="CS") {
      Ns=sum(nk*(nk-1))
    }else if (worktype=="AR"){
      Ns=sum(nk-1);Ns1=sum(nk-2)
    }else if (worktype=="Tri"){
      Ns=sum(nk-1);Ns1=sum(nk-2)
    }else{
      Ns=nsub
    }

    sumrd<-sumrd1<-0;R<-sumjk<-matrix(0,max(nk),max(nk));

    for (i in 1:nsub) {
      if (nk[i]>=2){
        za=sd_resid[(cn[i]+1):cn[i+1]]
        if (worktype=="CS"){
          for (j in 1:nk[i]){
            for (k in 1:nk[i]){
              if (k!=j) sumrd=sumrd+za[j]*za[k]
            }
          }
        }else if(worktype=="AR"){
          for (j in 2:nk[i]) {
            sumrd=sumrd+za[j]*za[j-1]
            if (j>=3) sumrd1=sumrd1+za[j]*za[j-2]
          }
        }else if(worktype=="Tri"){
          for (j in 2:nk[i]) sumrd=sumrd+za[j]*za[j-1]
        }else{
          for (j in 1:nk[i]){
            for (k in 1:nk[i]){
              R[j,k]=R[j,k]+za[j]*za[k]
              sumjk[j,k]=sumjk[j,k]+1
            }
          }
        }
      }
    }



    if (worktype=="CS"){
      af=sumrd/Ns;R=diag(max(nk))
      for (i in 1:max(nk)){
        for (j in 1:max(nk)){
          if (i!=j) R[i,j]=af
        }
      }
    }else if(worktype=="AR"){
      #af=(sumrd/Ns+sqrt(max(sumrd1,0))/Ns1)/2;R=diag(max(nk))
      af=sumrd/Ns;R=diag(max(nk))
      for (i in 1:max(nk)){
        for (j in 1:max(nk)){
          R[i,j]=af^(abs(i-j))
        }
      }
    }else if(worktype=="Tri"){
      af=(sumrd/Ns+sqrt(sumrd1)/Ns1)/2;
      R=diag(max(nk))
      for (i in 1:max(nk)){
        for (j in 1:max(nk)){
          if (abs(i-j)<=1) R[i,j]=af^(abs(i-j))
        }
      }
    }else if(worktype=="Ind"){
      R=diag(max(nk))
    }else{

      R=R/(nsub-nx);
      temp<-sqrt(diag(R)); R<-t(t(R/temp)/temp)

    }

    ############induce smoothing ##################
    U2=0;H2=0;
    r2=sqrt(diag(x%*%t(x)))
    mu2=x%*%beta;hh=y-mu2
    fr=(pnorm(nsub^0.5*(y-mu2)/r2,0,1)-rep(1-tau,N))


    for (i in 1:nsub){
      a_ind = (cn[i]+1):cn[i+1]
      Ra = R[1:nk[i],1:nk[i]];
      fra=fr[a_ind]
      r2a=r2[a_ind]
      ha=hh[a_ind]
      if (nk[i]==1){
        Ga=as.vector(f0[a_ind])
        xa = matrix(x[a_ind,],nrow=1);
        Ha2=as.vector(dnorm(nsub^0.5*ha/r2a,0,1)/r2a)
      }else{
        Ga=diag(as.vector(f0[a_ind]))
        xa = x[a_ind,];
        Ha2=diag(as.vector(dnorm(nsub^0.5*ha/r2a,0,1)/r2a))
      }
      U2=U2+t(xa)%*%Ga%*%ginv(Ra)%*%fra
      H2=H2+nsub^0.5*t(xa)%*%Ga%*%ginv(Ra)%*%Ha2%*%xa
    }



    if(lambda>0){
      #penalty
      eps <- 10^-6
      #scad
      pe <- as.vector(pp_scad(abs(as.vector(beta)),lambda)/(abs(as.vector(beta))+eps))
      if(length(pe)!=1){
        sigma<- nsub*diag(pe)
      }else{
        sigma <- nsub*pe
      }
      diff2=ginv(H2 + sigma)%*%(U2-sigma%*%beta)
    }else{
      diff2=ginv(H2)%*%U2
    }



    beta=beta+w*diff2
    ghat=x%*%beta
    mcl = mean(check_loss(y-ghat,tau))

    if (max(abs(diff2))>100) break

    diff2a=max(abs(w*diff2))
    w=w/2
  }


  converge = TRUE
  if (kk>=3 || iteration>=max_it || max(abs(diff2))>100) {
    #message("iterations: ",iteration)
    #message("max(abs(diff)): ",max(abs(diff2)))
    #message("current beta is",paste(round(beta,4),collapse = " "))
    converge = FALSE
  }

  if(!is.null(removed_ind)){
    beta_out = beta_all
    beta_out[removed_ind]=0
    beta_out[-removed_ind]=beta
    X_selected = setdiff((1:length(beta_all)),paste(removed_ind))
  }else{
    beta_out = beta_all
    X_selected = (1:length(beta_all))
  }

  #calculate hbic
  hbic = log(mcl*N)+(log(nsub)/(2*nsub))*log(log(NCOL(x)))*length(X_selected)
  list(beta=beta_out,g=ghat,R=R,X_selected=X_selected,mcl=mcl,hbic=hbic,converge=converge)

}

#Quantile Longitudinal GEE SCAD HBIC
qlingee_scad_hbic<-function(x,y,betaint,nk,worktype,f0=rep(1,length(y)),tau=0.5,lambda=NULL,max_it = 100,cutoff=10^-1){
  if(is.null(lambda)){
    #similiar to glmnet
    lambda_max=10
    lambda.min.ratio=ifelse(length(nk) > NCOL(x), 1e-03, 0.01 )
    lambda = exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio),
                     length.out = 30))
  }
  #library(doParallel)
  #n_core = parallel::detectCores(FALSE,FALSE)
  #cl <- parallel::makeCluster(n_core)
  #doParallel::registerDoParallel(cl)

  fit_results <- foreach::foreach(l = lambda,.combine=rbind,.packages = c("MASS"))%do%{
    source("./functions/plrq.R", local = TRUE)
    source("./functions/utils.R", local = TRUE)
    result <- qlingee_scad(x,y,betaint,nk,worktype,f0,tau,l,max_it,cutoff)
    c(l,result$mcl,result$hbic,paste0(result$X_selected,collapse = ""))
  }

  #parallel::stopCluster(cl)
  fit_results = as.data.frame(fit_results)
  fit_results[,3] = as.numeric(fit_results[,3])
  #print(fit_results)
  best_lambda <- lambda[which(fit_results[,3]==min(fit_results[,3]))]
  if(length(best_lambda)>1) {best_lambda=best_lambda[1]}
  qlingee_scad(x,y,betaint,nk,worktype,f0,tau,best_lambda,max_it,cutoff)
}

#Quantile Longitudinal GEE SCAD cross-validation
qlingee_scad_cv<-function(x,y,betaint,nk,worktype,f0=rep(1,length(y)),tau=0.5,lambda=NULL,max_it = 100,cutoff=10^-3,nfold=3){
  if(is.null(lambda)){
    #need to change
    lambda = seq(0.01,0.5,0.2)
  }
  for (l in lambda){
    #this method only applies to balanced data, nk are the same
    nobs = nk[1]
    nsub = length(nk)
    sampled_sind = suppressWarnings(split(sample(1:nsub),1:nfold))
    cl = c()
    for(fold in 1:nfold){
      test_sind = sampled_sind[[fold]]
      train_sind = setdiff(1:nsub,test_sind)
      test_ind = as.vector(sapply((test_sind-1)*nobs+1,function(x) x:(x+nobs-1)))
      train_ind = as.vector(sapply((train_sind-1)*nobs+1,function(x) x:(x+nobs-1)))

      result <- qlingee_scad(x[train_ind,],y[train_ind],betaint,rep(nobs,length(train_sind)),worktype,f0,tau,l,max_it,cutoff)
      cl = c(cl,check_loss(y[test_ind]-x[test_ind,]%*%result$beta,tau))
    }

    mcl = mean(cl)
    print(mcl)
  }

}

