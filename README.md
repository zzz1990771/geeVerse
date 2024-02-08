# geeVerse

`geeVerse` is an R package to provide computationally efficient implementations of penalized generalized estimating equations for any combination of 1) simultaneous variable selection and estimation for high and even ultra-high dimensional data, 2) conditional quantile or mean regression, and 3) longitudinal or cross-sectional data analysis.

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
qpgee.fit = qpgee(x,y,tau=0.5,nk=n_obs)
qpgee.fit$beta

#fit qpgee with auto selected lambda with parallel computing
qpgee.fit = qpgee_tune(x,y,tau=0.5,nk=n_obs,ncore=10)
qpgee.fit$beta
qpgee.fit$best_lambda
```

