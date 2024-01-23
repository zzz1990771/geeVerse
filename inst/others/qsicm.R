
Jbeta<-function(betaold){
betanew=c(sqrt(1-sum(betaold^2)),betaold)
J1=-betaold/betanew[1]
J=rbind(J1,diag(length(J1)))
}


quantsi<-function(x,y,z,thtaint,tau){

n=dim(x)[1];dx=dim(x)[2];dz=dim(z)[2]
thta=thtaint[-1]

iter<-0;diff=1;

while (diff>1e-3 & iter<40){

iter=iter+1;
betaold=thta[1:(dx-1)];af=thta[dx:(dx+dz-1)]
if (sum(betaold^2)>1) betaold=betaold*0.8;
betanew<-c(sqrt(1-sum(betaold^2)),betaold);
ul=x%*%betanew;B=NULL;
kl=min(ul)+(max(ul)-min(ul)+0.0001)*(0:(N+1))/(N+1)
kl=c(seq(min(ul),min(ul), length = (p-1)),kl,seq(max(ul),max(ul), length = (p-1)))#add knotes to both ends
B=splineDesign(knots=kl, x=ul,ord = p);
B1=splineDesign(knots=kl, x=ul,ord = p,derivs=rep(1,n))

delta<-coefficients(rq(y-z%*%af~0+B,tau=tau));
gz<-B%*%delta;gz1<-B1%*%delta
eta=ginv(t(B)%*%B)%*%t(B)%*%x;xt=x-B%*%eta
#etaz=ginv(t(B)%*%B)%*%t(B)%*%z;zt=z-B%*%etaz

thtaold=thta;
J=Jbeta(as.vector(betaold));
ya=y-gz+diag(as.vector(gz1))%*%xt%*%J%*%betaold;
xz=cbind(diag(as.vector(gz1))%*%xt%*%J,z)
thta=coefficients(rq(ya~0+xz,tau=tau))
diff=max(abs(thtaold-thta))
}

betaold=thta[1:(dx-1)];af=thta[dx:(dx+dz-1)]
if (sum(betaold^2)>1) betaold=betaold*0.8;
betanew<-c(sqrt(1-sum(betaold^2)),betaold)
return(list(beta=betanew,af=af))
}






qsicm_est<-function(x,y,z,beta,tau){
n=dim(x)[1];dx<-dim(x)[2];dz<-dim(z)[2]
betaint<-beta[-1]

###
############################################################################

    fn<-function(xx){
      #####
      betaold=xx
    	betaold1=sqrt(abs(1-sum(betaold^2)))
    	betaold=c(betaold1,betaold)
    	ul=x%*%betaold
    	B=NULL
    	kl=min(ul)+(max(ul)-min(ul)+0.0001)*(0:(N+1))/(N+1)#internal knot
      kl=c(seq(min(ul),min(ul), length = (p-1)),kl,seq(max(ul),max(ul), length = (p-1)))#add knotes to both ends
      B=splineDesign(knots=kl, x=ul,ord = p)

      Bz=NULL;
      for (j in 1:dz){
      Bz<-cbind(Bz,z[,j]*B)
      }

      lambdaold=coefficients(rq(y~Bz+0,tau=tau))
      res=y-(Bz%*%lambdaold)
      Ln=sum(res*(tau-(res<0)))
      Ln
    }
##############################################################################


##gr gradient of objective function
###############################################################################

gr <- function(xx){
      Ln1 <- rep(NA, length(xx))
      betaold=xx
    	betaold1=sqrt(abs(1-sum(betaold^2)))
    	betaold=c(betaold1,betaold)
    	ul=x%*%betaold
    	B=NULL

    	kl=min(ul)+(max(ul)-min(ul)+0.0001)*(0:(N+1))/(N+1)#internal knot
      kl=c(seq(min(ul),min(ul), length = (p-1)),kl,seq(max(ul),max(ul), length = (p-1)))#add knotes to both ends
      B=splineDesign(knots=kl, x=ul,ord = p)

      Bz=NULL;
      for (j in 1:dz){
      Bz<-cbind(Bz,z[,j]*B)
      }
      lambdaold=coefficients(rq(y~Bz+0,tau=tau))
      res=y-(Bz%*%lambdaold)

      J=Jbeta(xx)
      eta=ginv(t(Bz)%*%Bz)%*%t(Bz)%*%x
      xtilda=x-Bz%*%eta

      kl=min(ul)+(max(ul)-min(ul)+0.0001)*(0:(N+1))/(N+1)
      kl=c(seq(min(ul),min(ul), length = (p-1)),kl,seq(max(ul),max(ul), length = (p-1)))
      Blprime=splineDesign(knots=kl, x=ul,ord = p,derivs=rep(1,n))
      Blprimez=NULL;
      for (j in 1:dz){
        Blprimez<-cbind(Blprimez,z[,j]*Blprime)
      }
      mhatlprime=Blprimez%*%lambdaold

    Ln1=t(J)%*%t(xtilda)%*%diag(as.vector(mhatlprime))%*%as.vector(tau-(res<0))
    #Ln1=t(J)%*%t(x)%*%diag(as.vector(mhatlprime))%*%as.vector(res)

    Ln1=-Ln1
    Ln1
   }

########################################################################################

      #constraint matrix
      hin <- function(xx) {
      h <- rep(NA, 1)
      h[1]<- (1-sum(xx^2))
      }


     hin.jac <- function(xx) {
          j <- matrix(NA, 1, length(xx))
          j[1,]<--2*xx
          j
      }

    thtanew=constrOptim.nl(par=betaint, fn=fn,gr=gr,hin=hin, hin.jac=hin.jac)

    betanew=thtanew$par
    betanew1=sqrt(abs(1-sum(betanew^2)))
    betanew=c(betanew1,betanew)
    return(list(betanew=betanew))
}


lssicm_est<-function(x,y,z,beta){
n=dim(x)[1];dx<-dim(x)[2];dz<-dim(z)[2]
betaint<-beta[-1]

###fn use ls
############################################################################

    fn<-function(xx){
      #####remove one and add back
      betaold=xx
    	betaold1=sqrt(abs(1-sum(betaold^2)))
    	betaold=c(betaold1,betaold)
    	ul=x%*%betaold
    	B=NULL
    	kl=min(ul)+(max(ul)-min(ul)+0.0001)*(0:(N+1))/(N+1)#internal knot
      kl=c(seq(min(ul),min(ul), length = (p-1)),kl,seq(max(ul),max(ul), length = (p-1)))#add knotes to both ends
      B=splineDesign(knots=kl, x=ul,ord = p)

      Bz=NULL;
      for (j in 1:dz){
      Bz<-cbind(Bz,z[,j]*B)
      }

      lambdaold=coefficients(lm(y~Bz+0))
      res=y-(Bz%*%lambdaold)
      Ln=0.5*sum(res^2)
      Ln
    }
##############################################################################


##gr gradient of objective function
###############################################################################

gr <- function(xx){
      Ln1 <- rep(NA, length(xx))
      betaold=xx
    	betaold1=sqrt(abs(1-sum(betaold^2)))
    	betaold=c(betaold1,betaold)
    	ul=x%*%betaold
    	B=NULL

    	kl=min(ul)+(max(ul)-min(ul)+0.0001)*(0:(N+1))/(N+1)#internal knot
      kl=c(seq(min(ul),min(ul), length = (p-1)),kl,seq(max(ul),max(ul), length = (p-1)))#add knotes to both ends
      B=splineDesign(knots=kl, x=ul,ord = p)

      Bz=NULL;
      for (j in 1:dz){
      Bz<-cbind(Bz,z[,j]*B)
      }
      lambdaold=coefficients(lm(y~Bz+0))
      res=y-(Bz%*%lambdaold)

      J=Jbeta(xx)
      eta=ginv(t(Bz)%*%Bz)%*%t(Bz)%*%x
      xtilda=x-Bz%*%eta

      kl=min(ul)+(max(ul)-min(ul)+0.0001)*(0:(N+1))/(N+1)
      kl=c(seq(min(ul),min(ul), length = (p-1)),kl,seq(max(ul),max(ul), length = (p-1)))
      Blprime=splineDesign(knots=kl, x=ul,ord = p,derivs=rep(1,n))
      Blprimez=NULL;
      for (j in 1:dz){
        Blprimez<-cbind(Blprimez,z[,j]*Blprime)
      }
      mhatlprime=Blprimez%*%lambdaold

    #Ln1=t(J)%*%t(xtilda)%*%diag(as.vector(mhatlprime))%*%as.vector(tau-(res<0))
    Ln1=t(J)%*%t(xtilda)%*%diag(as.vector(mhatlprime))%*%as.vector(res)

    Ln1=-Ln1
    Ln1
   }

########################################################################################

      #constraint matrix
      hin <- function(xx) {
      h <- rep(NA, 1)
      h[1]<- (1-sum(xx^2))
      }


     hin.jac <- function(xx) {
          j <- matrix(NA, 1, length(xx))
          j[1,]<--2*xx
          j
      }

    thtanew=constrOptim.nl(par=betaint, fn=fn,gr=gr,hin=hin, hin.jac=hin.jac)

    betanew=thtanew$par
    betanew1=sqrt(abs(1-sum(betanew^2)))
    betanew=c(betanew1,betanew)
    return(list(betanew=betanew))
}



quantile_est<-function(x,y,beta,tau){
n=dim(x)[1];dx<-dim(x)[2]
betaint<-beta[-1]
N=2;p=4;

###fn use ls
############################################################################

    fn<-function(xx){
      #####remove one and add back
      betaold=xx
    	betaold1=sqrt(abs(1-sum(betaold^2)))
    	betaold=c(betaold1,betaold)
    	ul=x%*%betaold
    	B=NULL

    	kl=min(ul)+(max(ul)-min(ul)+0.0001)*(0:(N+1))/(N+1)#internal knot
      kl=c(seq(min(ul),min(ul), length = (p-1)),kl,seq(max(ul),max(ul), length = (p-1)))#add knotes to both ends
      B=splineDesign(knots=kl, x=ul,ord = p)

      lambdaold=coefficients(rq(y~B+0,tau=tau))
      res=y-(B%*%lambdaold)
      Ln=sum(res*(tau-(res<0)))
      Ln
    }
##############################################################################


##gr gradient of objective function
###############################################################################

gr <- function(xx){
      Ln1 <- rep(NA, length(xx))
      betaold=xx
    	betaold1=sqrt(abs(1-sum(betaold^2)))
    	betaold=c(betaold1,betaold)
    	ul=x%*%betaold
    	B=NULL

    	kl=min(ul)+(max(ul)-min(ul)+0.0001)*(0:(N+1))/(N+1)#internal knot
      kl=c(seq(min(ul),min(ul), length = (p-1)),kl,seq(max(ul),max(ul), length = (p-1)))#add knotes to both ends
      B=splineDesign(knots=kl, x=ul,ord = p)

      lambdaold=coefficients(rq(y~B+0,tau=tau))
      res=y-(B%*%lambdaold)

     J=Jbeta(xx)
     eta=ginv(t(B)%*%B)%*%t(B)%*%x
     xtilda=x-B%*%eta

    kl=min(ul)+(max(ul)-min(ul)+0.0001)*(0:(N+1))/(N+1)
    kl=c(seq(min(ul),min(ul), length = (p-1)),kl,seq(max(ul),max(ul), length = (p-1)))
    Blprime=splineDesign(knots=kl, x=ul,ord = p,derivs=rep(1,n))
    mhatlprime=Blprime%*%lambdaold

    Ln1=t(J)%*%t(xtilda)%*%diag(as.vector(mhatlprime))%*%as.vector(tau-(res<0))
    #Ln1=t(J)%*%t(x)%*%diag(as.vector(mhatlprime))%*%as.vector(res)

    Ln1=-Ln1
    #Ln1
   }

########################################################################################

      #constraint matrix
      hin <- function(xx) {
      h <- rep(NA, 1)
      h[1]<- (1-sum(xx^2))
      }


     hin.jac <- function(xx) {
          j <- matrix(NA, 1, length(xx))
          j[1,]<--2*xx
          j
      }

    thtanew=constrOptim.nl(par=betaint, fn=fn,gr=gr,hin=hin, hin.jac=hin.jac)

    betanew=thtanew$par
    betanew1=sqrt(abs(1-sum(betanew^2)))
    betanew=c(betanew1,betanew)
    return(list(betanew=betanew))
}


qplsim_est<-function(x,z,y,thta,tau){
n<-dim(x)[1];dx<-dim(x)[2];dz<-dim(z)[2]
beta=thta[1:dx];af=thta[(dx+1):(dx+dz)];
betaint<-beta[-1];thtaint=c(betaint,af)

###fn use ls
############################################################################

    fn<-function(xx){
      #####remove one and add back
      betaold=xx[1:(dx-1)];af=xx[dx:(dx+dz-1)]
    	betaold1=sqrt(abs(1-sum(betaold^2)))
    	betaold=c(betaold1,betaold)
    	ul=x%*%betaold
    	B=NULL

    	kl=min(ul)+(max(ul)-min(ul)+0.0001)*(0:(N+1))/(N+1)#internal knot
      kl=c(seq(min(ul),min(ul), length = (p-1)),kl,seq(max(ul),max(ul), length = (p-1)))#add knotes to both ends
      B=splineDesign(knots=kl, x=ul,ord = p)

      lambdaold=coefficients(rq(y-z%*%af~B+0,tau=tau))
      res=y-(B%*%lambdaold)-z%*%af
      Ln=sum(res*(tau-(res<0)))
      Ln
    }
##############################################################################


##gr gradient of objective function
###############################################################################

gr <- function(xx){
      Ln1 <- rep(NA, length(xx))
      betaold=xx[1:(dx-1)];af=xx[dx:(dx+dz-1)]
    	betaold1=sqrt(abs(1-sum(betaold^2)))
    	betaold=c(betaold1,betaold)
    	ul=x%*%betaold
    	B=NULL

    	kl=min(ul)+(max(ul)-min(ul)+0.0001)*(0:(N+1))/(N+1)#internal knot
      kl=c(seq(min(ul),min(ul), length = (p-1)),kl,seq(max(ul),max(ul), length = (p-1)))#add knotes to both ends
      B=splineDesign(knots=kl, x=ul,ord = p)

      lambdaold=coefficients(rq(y-z%*%af~B+0,tau=tau))
      res=y-(B%*%lambdaold)-z%*%af
      J=Jbeta(xx[1:(dx-1)])
      eta=ginv(t(B)%*%B)%*%t(B)%*%x
      etaz=ginv(t(B)%*%B)%*%t(B)%*%x
      zt=z-B%*%etaz

    kl=min(ul)+(max(ul)-min(ul)+0.0001)*(0:(N+1))/(N+1)
    kl=c(seq(min(ul),min(ul), length = (p-1)),kl,seq(max(ul),max(ul), length = (p-1)))
    Blprime=splineDesign(knots=kl, x=ul,ord = p,derivs=rep(1,n))
    mhatlprime=Blprime%*%lambdaold

    Ln1=rbind(t(J)%*%t(xtilda)%*%diag(as.vector(mhatlprime)),t(zt))%*%as.vector(tau-(res<0))
    #Ln1=t(J)%*%t(x)%*%diag(as.vector(mhatlprime))%*%as.vector(res)

    Ln1=-Ln1
    Ln1
   }

########################################################################################

      #constraint matrix
      hin <- function(xx) {
      h <- rep(NA, 1)
      h[1]<- (1-sum(xx^2))
      }


     hin.jac <- function(xx) {
          j <- matrix(NA, 1, length(xx))
          j[1,]<--2*xx
          j
      }

    thtanew=constrOptim.nl(par=thtaint, fn=fn,gr=gr,hin=hin, hin.jac=hin.jac)

    thtanew=thtanew$par;betanew1=thtanew[1:(dx-1)]
    betanew1=sqrt(abs(1-sum(betanew^2)))
    betanew=c(betanew1,betanew)
    return(list(betanew=betanew,af=thtanew[(dx-1):(dx+dz-1)]))
}



ls_est<-function(x,y,beta){
n=dim(x)[1];dx<-dim(x)[2]
betaint<-beta[-1]
N=2;p=4;

###fn use ls
############################################################################

    fn<-function(xx){
      #####remove one and add back
      betaold=xx
    	betaold1=sqrt(abs(1-sum(betaold^2)))
    	betaold=c(betaold1,betaold)
    	ul=x%*%betaold
    	B=NULL

    	kl=min(ul)+(max(ul)-min(ul)+0.0001)*(0:(N+1))/(N+1)#internal knot
      kl=c(seq(min(ul),min(ul), length = (p-1)),kl,seq(max(ul),max(ul), length = (p-1)))#add knotes to both ends
      B=splineDesign(knots=kl, x=ul,ord = p)

      lambdaold=coefficients(lm(y~B+0))
      res=y-(B%*%lambdaold)
      Ln=0.5*sum(res^2)
      Ln
    }
##############################################################################


##gr gradient of objective function
###############################################################################

gr <- function(xx){
      Ln1 <- rep(NA, length(xx))
      betaold=xx
    	betaold1=sqrt(abs(1-sum(betaold^2)))
    	betaold=c(betaold1,betaold)
    	ul=x%*%betaold
    	B=NULL

    	kl=min(ul)+(max(ul)-min(ul)+0.0001)*(0:(N+1))/(N+1)#internal knot
      kl=c(seq(min(ul),min(ul), length = (p-1)),kl,seq(max(ul),max(ul), length = (p-1)))#add knotes to both ends
      B=splineDesign(knots=kl, x=ul,ord = p)

      lambdaold=coefficients(lm(y~B+0))
      res=y-(B%*%lambdaold)

     J=Jbeta(xx)
     eta=ginv(t(B)%*%B)%*%t(B)%*%x
     xtilda=x-B%*%eta

    kl=min(ul)+(max(ul)-min(ul)+0.0001)*(0:(N+1))/(N+1)
    kl=c(seq(min(ul),min(ul), length = (p-1)),kl,seq(max(ul),max(ul), length = (p-1)))
    Blprime=splineDesign(knots=kl, x=ul,ord = p,derivs=rep(1,n))
    mhatlprime=Blprime%*%lambdaold

    #Ln1=t(J)%*%t(xtilda)%*%diag(as.vector(mhatlprime))%*%as.vector(tau-(res<0))
    Ln1=t(J)%*%t(xtilda)%*%diag(as.vector(mhatlprime))%*%as.vector(res)

    Ln1=-Ln1
    Ln1
   }

########################################################################################

      #constraint matrix
      hin <- function(xx) {
      h <- rep(NA, 1)
      h[1]<- (1-sum(xx^2))
      }


     hin.jac <- function(xx) {
          j <- matrix(NA, 1, length(xx))
          j[1,]<--2*xx
          j
      }

    thtanew=constrOptim.nl(par=betaint, fn=fn,gr=gr,hin=hin, hin.jac=hin.jac)

    betanew=thtanew$par
    betanew1=sqrt(abs(1-sum(betanew^2)))
    betanew=c(betanew1,betanew)
    return(list(betanew=betanew))
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

    eta=geninv(t(B)%*%B)%*%t(B)%*%x;xb=x-B%*%eta
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
