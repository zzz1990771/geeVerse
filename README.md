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
sim_data <- generateData(nsub = 20, nobs = rep(10, 20),  p = 20,
                         beta0 = c(rep(1,5),rep(0,15)), rho = 0.1, correlation = "AR1",
                          dis = "normal", ka = 1)

X=sim_data$X
y=sim_data$y

#fit qpgee with auto selected lambda
qpgee.fit = qpgee(X,y,tau=0.5,nobs=rep(10, 20),ncore=1)
qpgee.fit$beta
```

