qlinqif_scad<-function(x,y,betaint,nk,worktype,f0=rep(1,length(y)),tau=0.5,lambda=0.1,max_it = 100,cutoff=10^-1){
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

    ############induce smoothing ##################
    U2=0;H2=0;
    U2C=0;
    r2=sqrt(diag(x%*%t(x)))
    mu2=x%*%beta;hh=y-mu2
    fr=(pnorm(nsub^0.5*(y-mu2)/r2,0,1)-rep(1-tau,N))


    for (i in 1:nsub){
      #index for x
      a_ind = (cn[i]+1):cn[i+1]
      fra=fr[a_ind]
      r2a=r2[a_ind]
      ha=hh[a_ind]

      corstr = worktype
      #qif basis
      ni <- nk[i]
      m0 <- diag(ni)
      if (corstr == "Ind") {
        m1 <- matrix(rep(0,ni*ni),ni)
      }else if (corstr == "CS") {
        m1 <- matrix(rep(1,ni*ni),ni) - m0
      }else if (corstr == "AR") {
        m1 <- matrix(rep(0,ni*ni),ni)
        for (k in 1:ni) {
          for (l in 1:ni) {
            if (abs(k-l)==1) m1[k,l] <-1
          }
        }
      }


      #set
      if (nk[i]==1){
        Ga=as.vector(f0[a_ind])
        xa = matrix(x[a_ind,],nrow=1);
        Ha2=as.vector(dnorm(nsub^0.5*ha/r2a,0,1)/r2a)
      }else{
        Ga=diag(as.vector(f0[a_ind]))
        xa = x[a_ind,];
        Ha2=diag(as.vector(dnorm(nsub^0.5*ha/r2a,0,1)/r2a))
      }

      #
      if (corstr == "Ind") {
        s_m0=t(xa)%*%Ga%*%m0%*%fra
        U2=U2 + s_m0
        U2C = U2C + s_m0 %*% t(s_m0)
        h_m0 = nsub^0.5*t(xa)%*%Ga%*%m0%*%Ha2%*%xa
        H2 = H2 + h_m0
      }else {
        s_m0=t(xa)%*%Ga%*%m0%*%fra
        s_m1=t(xa)%*%Ga%*%m1%*%fra
        U2 = U2 + c(s_m0,s_m1)
        U2C = U2C + c(s_m0,s_m1) %*% t(c(s_m0,s_m1))

        h_m0 = nsub^0.5*t(xa)%*%Ga%*%m0%*%Ha2%*%xa
        h_m1 = nsub^0.5*t(xa)%*%Ga%*%m1%*%Ha2%*%xa
        H2 = H2 + rbind(h_m0,h_m1)

      }


    }


    #scale
    #H2 = H2/nsub
    #U2 = U2/nsub
    #U2C = U2C/nsub

    U2 = as.matrix(U2,nrow=length(U2))
    #H2 arsumgfirstdev, U2 arsumg U2C arsumc
    arqif1dev <- t(H2) %*% ginv(U2C) %*% U2
    arqif2dev <- t(H2) %*% ginv(U2C) %*% H2


    if(lambda>0){
      #penalty
      eps <- 10^-6
      #scad
      pe <- as.vector(pp_scad_sim(abs(as.vector(beta)),lambda)/(abs(as.vector(beta))+eps))
      if(length(pe)!=1){
        sigma<- nsub*diag(pe)
      }else{
        sigma <- nsub*pe
      }
      diff2=ginv(arqif2dev + sigma)%*%(arqif1dev-sigma%*%beta)
    }else{
      diff2=ginv(arqif2dev)%*%arqif1dev
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


qlinqif_scad_hbic<-function(x,y,betaint,nk,worktype,f0=rep(1,length(y)),tau=0.5,lambda=NULL,max_it = 100,cutoff=10^-1){
  if(is.null(lambda)){
    #similiar to glmnet
    lambda_max=1
    lambda.min.ratio=ifelse(length(nk) > NCOL(x), 1e-03, 0.01 )
    lambda = exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio),
                     length.out = 30))
  }
  library(doParallel)
  n_core = parallel::detectCores(FALSE,FALSE)
  cl <- parallel::makeCluster(n_core)
  doParallel::registerDoParallel(cl)

  fit_results <- foreach::foreach(l = lambda,.combine=rbind,.packages = c("MASS"))%dopar%{
    source("./R/plrq.R", local = TRUE)
    source("./R/qqif.R", local = TRUE)
    source("./R/utils.R", local = TRUE)
    result <- qlinqif_scad(x,y,betaint,nk,worktype,f0,tau,l,max_it,cutoff)
    c(l,result$mcl,result$hbic,paste0(result$X_selected,collapse = ""))
  }

  parallel::stopCluster(cl)
  fit_results = as.data.frame(fit_results)
  fit_results[,3] = as.numeric(fit_results[,3])
  print(fit_results)
  best_lambda <- lambda[which(fit_results[,3]==min(fit_results[,3]))]
  if(length(best_lambda)>1) {best_lambda=best_lambda[1]}
  qlinqif_scad(x,y,betaint,nk,worktype,f0,tau,best_lambda,max_it,cutoff)
}
