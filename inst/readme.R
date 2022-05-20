# README / TODO

# This file describes how to generate the numerical results in the paper.

#### Autoregressive process of order 1

# The files starting with 'ar1' refer to the autoregressive process of order 1,
# for which the asymptotic variance is known in closed form.

# ar1functions.R: defines the model and algorithms
# ar1asymptvarestimation.R: run to generate ar1.uavar.RData.
source("inst/ar1asymptvarestimation.R")

# ar1fishyestimation.R: tun to generate ar1.fishyfunction.RData
source("inst/ar1fishyestimation.R")

# ar1longrun.R: creates ar1.longruns.RData (large file)
# source("inst/ar1longrun.R")

# ar1batchmeans.R: creates ar1batchmeans.RData
# source("inst/ar1batchmeans.R")
# ar1spectralvar.R: creates ar1spectralvar.RData
# source("inst/ar1spectralvar.R")

# ar1plots.R: creates figures and tables
source("inst/ar1plots.R")

# ar1boundedh.R: considers a bounded test function
# ar1unbiasedmcmc.R
# other files in backup/

#### Cauchy-Normal model

# The files starting with 'cauchynormal' refer to the posterior in a
# Cauchy-Normal model used in a paper by Christian Robert,
# "Convergence control methods for Markov chain Monte Carlo algorithms", 1995.

# cauchynormalfunctions.R: defines the model and two MCMC algorithms targeting the posterior

# cauchynormalfishyestimation.R: to generate cauchynormal.fishyfunction.RData
source("inst/cauchynormalfishyestimation.R")

# cauchynormalplots.R: creates all plots and tables for this model
source("inst/cauchynormalplots.R")

## other file: cauchynormallongrun.R, cauchynormalunbiasedmcmc.R

#### High-dimensional Bayesian Regression
# The files starting with 'highdimreg' refer to an example of 
# Bayesian linear regression with a shrinkage prior, 
# the coupling of which is the topic of the article
# "Coupled Markov chain Monte Carlo for high-dimensional regression with Half-t priors"
# Niloy Biswas, Anirban Bhattacharya, Pierre E Jacob, James E Johndrow, 2022.

# highdimregfunctions.R: defines the Gibbs sampler and its coupling

# highdimregavarestimation.R: generates highdimreg.meetings.RData and highdimreg.uavar.RData
source("inst/highdimregavarestimation.R")

# highdimreglongrun.R: generates highdimreg.longrun.RData
source("inst/highdimreglongrun.R")

# highdimregplots.R: generates tables and plots
source("inst/highdimregplots.R")

# highdimreg_test.R, highdimreg_avarestimation_mcmc.R,
# highdimregunbiasedmcmc.R



