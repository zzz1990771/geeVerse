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
sim_data <- generate_data(
  nsub = 50, nobs = rep(5, 50), p = 10,
  beta0 = c(rep(1, 5), rep(0, 5)), rho = 0.3
)

# 2. Fit the model using the formula interface
fit <- qpgee(
  y ~ . - id,
  data = sim_data,
  id = sim_data$id,
  tau = 0.5,
  method = "HBIC"
)

# 3. View the summary of the results
summary(fit)
```

## Updates
This package was re-factored with major functions to make it more consistent with other R packages.
