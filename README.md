
## unbiasedpoisson

This package contains scripts that generate the figures and tables of
the article “Solving the Poisson equation using coupled Markov chains”
by Randal Douc, Pierre E. Jacob, Anthony Lee and Dootika Vats.

### Installation

The package can be installed from R via:

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("pierrejacob/unbiasedpoisson")
```

This should install `Rcpp`, `RcppEigen`, `tictoc` automatically if
needed.

Most scripts depend on packages which can be installed via:

``` r
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("doParallel")) install.packages("doParallel")
if (!require("doRNG")) install.packages("doRNG")
if (!require("latex2exp")) install.packages("latex2exp")
if (!require("ggthemes")) install.packages("ggthemes")
```

Some scripts also depend on the following packages:

``` r
if (!require("boot")) install.packages("boot")
if (!require("mcmcse")) install.packages("mcmcse")
if (!require("fftwtools")) install.packages("fftwtools")
if (!require("RcppArmadillo")) install.packages("RcppArmadillo")
```
