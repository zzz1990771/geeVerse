# geeVerse

`geeVerse` is an R package to provide computationally efficient implementations of penalized generalized estimating equations for conditional quantiles and the conditional mean, when the risk factors are high dimensional and the data is longitudinal.

## Installation

You can install the latest version of `geeVerse` from GitHub with:

```R
# install.packages("devtools")
devtools::install_github("zzz1990771/geeVerse")
```

## Usage and Example:

After installation, you can load the package as usual:

```R
library(geeVerse)
```

To get detailed documentation on the `qpgee` function, use:

```R
?qpgee
```

This will show you the function's usage, arguments, and examples.

Running an Example:

```R
#settings
n=n_sub=400
p=200
beta0=rep(1,7)
p0=length(beta0)
beta = c(beta0,rep(0,p-p0))
n_obs<-rep(10,n_sub);
ka = 1
rho=0.6
type="ar"
dis="normal"
n_sim = 2

#generate errors for each subject
e = NULL
id<-NULL
for (i in 1:n_sub){
  id<-c(id,rep(i,n_obs[i]))
  sigmai=Siga_cov(rho,type,n_obs[i])
  if (dis=="normal") ei=MASS::rmvnorm(1, mean=rep(0, n_obs[i]), sigma=sigmai)
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
qpgee(x,y,tau=0.5,nk=n_obs)
```

