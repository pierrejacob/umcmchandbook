
## umcmchandbook

This package contains scripts that generate the figures of the document
entitled “Unbiased Markov Chain Monte Carlo: what, why and how”, by Yves
F. Atchadé and Pierre E. Jacob, 2023.

See companion website at <https://pierrejacob.quarto.pub/unbiased-mcmc>.

### Installation

The package can be installed from R via:

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("pierrejacob/umcmchandbook")
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
