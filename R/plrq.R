#' Quantile Penalized Generalized Estimating Equations (QPGEE)
#'
#' This function implements Quantile Penalized Generalized Estimating Equations
#' (QPGEE) for longitudinal data analysis. It estimates parameters using a
#' penalized quantile regression approach within a GEE framework, allowing for
#' different working correlation structures.
#'
#' @param x A matrix of predictors.
#' @param y A numeric vector of response variables.
#' @param tau The quantile to be estimated (default is 0.5, the median).
#' @param nk A numeric vector indicating the number of observations per subject.
#' @param worktype A string specifying the working correlation structure.
#'        Options include "CS" (Compound Symmetry), "AR" (Autoregressive),
#'        "Tri" (Tri-diagonal), and "Ind" (Independent).
#' @param lambda The penalty parameter for regularization (default is 0.1).
#' @param betaint Initial values for the beta coefficients. If NULL,
#'        non-longitudinal quantile regression is used for initialization.
#' @param max_it Maximum number of iterations (default is 100).
#' @param cutoff Threshold for coefficient shrinkage (default is 0.1).
#' @return A list containing the following components:
#'         \itemize{
#'           \item{beta}{Estimated beta coefficients.}
#'           \item{g}{Fitted values of the linear predictor.}
#'           \item{R}{Estimated working correlation matrix.}
#'           \item{X_selected}{Indices of selected predictors.}
#'           \item{mcl}{Mean check loss.}
#'           \item{hbic}{Hannan-Quinn Information Criterion value.}
#'           \item{converge}{Boolean indicating whether the algorithm converged.}
#'         }
#' @examples
#' # Example usage:
#'
#' #data generation settings
#' n=n_sub=400
#' p=200
#' beta0=rep(1,7)
#' p0=length(beta0)
#' beta = c(beta0,rep(0,p-p0))
#' n_obs<-rep(10,n_sub);
#' ka = 1
#' rho=0.6
#' type="ar"
#' dis="normal"
#'
#' #generate errors for each subject
#' e = NULL
#' id<-NULL
#' for (i in 1:n_sub){
#'   id<-c(id,rep(i,n_obs[i]))
#'   sigmai=Siga_cov(rho,type,n_obs[i])
#'   if (dis=="normal") ei=MASS::rmvnorm(1, mean=rep(0, n_obs[i]), sigma=sigmai)
#'   if (dis=="t") ei=mvtnorm::rmvt(1, sigmai, df = 4, delta = rep(0, n_obs[i]))
#'   e=c(e,ei);
#' }
#'
#' #generate y and X
#' N=sum(n_obs)
#' nk=n_obs
#' cn = c(0, cumsum(n_obs))
#' x=X=matrix(rnorm(N*p),N,p)
#' y=X%*%beta+(1+ka*abs(X[,1]))*e
#'
#' #fit qpgee
#' qpgee(x,y,tau=0.5,nk=n_obs)
#'
#' @export
#' @importFrom stats rq
#' @importFrom MASS ginv
#' @importFrom stats pnorm
#' @importFrom stats dnorm
#' @importFrom base diag
#' @importFrom base matrix
#' @importFrom base setdiff
#' @importFrom base abs
#' @importFrom base sum
#' @importFrom base rep
#' @importFrom base which
#' @importFrom base is.null
#' @importFrom base c
#' @importFrom base length
#' @importFrom base sum
#' @importFrom base cumsum
#' @importFrom base coefficients
#' @importFrom base max
#' @importFrom base t
#' @importFrom base sqrt
#' @importFrom base log
#' @importFrom base ncol
qpgee<-function(x,y,tau=0.5,nk=rep(1,length(y)),worktype="CS",lambda=0.1,betaint=NULL,f0=NULL,max_it = 100,cutoff=10^-1){
  x_all = x
  cn = c(0, cumsum(nk));nsub=length(nk);nx=dim(x)[2];
  N=sum(nk);


  #if initial beta is not provided, use non-longitudinal quantile regression as initial
  if(is.null(betaint)){
    betaint = coefficients(rq(y~0+X,tau = tau))
  }
  betaold=beta_all=beta=betaint;

  if(is.null(f0)){
    f0=rep(1,length(y))
  }


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
      pe <- as.vector(pp_scad_sim(abs(as.vector(beta)),lambda)/(abs(as.vector(beta))+eps))
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


qlingee_scad_hbic<-function(x,y,betaint,nk,worktype,f0=rep(1,length(y)),tau=0.5,lambda=NULL,max_it = 100,cutoff=10^-1){
  if(is.null(lambda)){
    #similiar to glmnet
    lambda_max=10
    lambda.min.ratio=ifelse(length(nk) > NCOL(x), 1e-03, 0.01 )
    lambda = exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio),
                     length.out = 30))
  }
  library(doParallel)
  n_core = parallel::detectCores(FALSE,FALSE)
  cl <- parallel::makeCluster(n_core)
  doParallel::registerDoParallel(cl)
  print(lambda)
  fit_results <- foreach::foreach(l = lambda,.combine=rbind,.packages = c("MASS"))%dopar%{
    source("../R/plrq.R", local = TRUE)
    source("../R/utils.R", local = TRUE)
    result <- qlingee_scad(x,y,betaint,nk,worktype,f0,tau,l,max_it,cutoff)
    c(l,result$mcl,result$hbic,paste0(result$X_selected,collapse = ""))
  }

  parallel::stopCluster(cl)
  fit_results = as.data.frame(fit_results)
  fit_results[,3] = as.numeric(fit_results[,3])
  print(fit_results)
  best_lambda <- lambda[which(fit_results[,3]==min(fit_results[,3]))]
  print(best_lambda)
  if(length(best_lambda)>1) {best_lambda=best_lambda[1]}
  qlingee_scad(x,y,betaint,nk,worktype,f0,tau,best_lambda,max_it,cutoff)
}

qlingee_scad_cv<-function(x,y,betaint,nk,worktype,f0,tau=0.5,lambda=NULL,max_it = 100,cutoff=10^-3,nfold=3){
  if(is.null(lambda)){
    #need to change
    lambda = seq(0.01,0.5,0.2)
  }
  for (l in lambda){
    #this method only applies to balanced data, nk are the same
    n_obs = nk[1]
    n_sub = length(nk)
    sampled_sind = suppressWarnings(split(sample(1:n_sub),1:nfold))
    cl = c()
    for(fold in 1:nfold){
      test_sind = sampled_sind[[fold]]
      train_sind = setdiff(1:n_sub,test_sind)
      test_ind = as.vector(sapply((test_sind-1)*n_obs+1,function(x) x:(x+n_obs-1)))
      train_ind = as.vector(sapply((train_sind-1)*n_obs+1,function(x) x:(x+n_obs-1)))

      result <- qlingee_scad(x[train_ind,],y[train_ind],betaint,rep(n_obs,length(train_sind)),worktype,f0,tau,l,max_it,cutoff)
      cl = c(cl,check_loss(y[test_ind]-x[test_ind,]%*%result$beta,tau))
    }

    mcl = mean(cl)
    print(mcl)
  }

}


lssimgee<-function(x,y,betaint,nk,id,worktype="exchangeable"){

  cn = c(0, cumsum(nk));
  nsub=length(nk);nx=dim(x)[2];
  N=sum(nk);
  beta=betaint
  betaold=beta[-1]


  K=2;p=4;J=K+p;B<-matrix(0,N,J)
  ul=x%*%beta;kl=min(ul)+(max(ul)-min(ul)+0.0001)*(0:(K+1))/(K+1)
  kl=c(seq(min(ul),min(ul), length = (p-1)),kl,seq(max(ul),max(ul), length = (p-1)))
  B=splineDesign(knots=kl, x=ul,ord = p);
  B1=splineDesign(knots=kl, x=ul,ord = p,derivs=rep(1,N))

  thta=coefficients(lm(y~0+B))
  ghat=B%*%thta;ghat1=B1%*%thta;




  betadiff = 1; iteration = 0;w = 1;kk=0
  diff1a=diff2a=1


  while (max(abs(diff1a))>1e-4 & max(abs(diff2a))>0.01 & iteration<100){

    iteration=iteration+1
    Jb=Jbeta(as.vector(betaold))

    eta=ginv(t(B)%*%B)%*%t(B)%*%x;xb=x-B%*%eta
    xs=diag(as.vector(ghat1))%*%xb%*%Jb
    ys=y-ghat+xs%*%betaold

    aa1<-gee(ys ~ 0+xs, id=id, corstr=worktype)
    #aa1<-gee(ys ~ 0+xs, id=id, corstr="AR-M",Mv=1)
    betaold0=betaold
    betaold=aa1$coefficients
    beta=c(sqrt(1-sum(betaold^2)),betaold)


    K=2;p=4;J=K+p;B<-matrix(0,N,J)
    ul=x%*%beta;kl=min(ul)+(max(ul)-min(ul)+0.0001)*(0:(K+1))/(K+1)
    kl=c(seq(min(ul),min(ul), length = (p-1)),kl,seq(max(ul),max(ul), length = (p-1)))
    B=splineDesign(knots=kl, x=ul,ord = p);
    B1=splineDesign(knots=kl, x=ul,ord = p,derivs=rep(1,N))

    #aa2<-gee(ys ~ 0+xs, id=id, corstr="exchangeable")
    aa2<-gee(y ~ 0+B, id=id, corstr=worktype)#corstr="AR-M",Mv=1)

    thtaold=thta;thta<-aa2$coefficients
    ghat=B%*%thta;ghat1=B1%*%thta;

    diff1a=max(abs(betaold0-betaold))
    diff2a=max(abs(thta-thtaold))

  }


  list(beta=beta,thta=thta,g=ghat)

}

#need further work on this
predict.plrq <- function (object, newdata)
  {
    if (missing(newdata)){
      return(napredict(object$na.action, object$fitted))
    }else {
      X <- x
    }
    pred <- drop(X %*% object$beta)
    pred
  }
