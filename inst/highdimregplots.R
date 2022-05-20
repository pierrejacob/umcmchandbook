rm(list = ls())
library(unbiasedpoisson)
setmytheme()
library(doParallel)
library(doRNG)
library(dplyr)
# register parallel cores
registerDoParallel(cores = 6)
# set RNG seed
set.seed(1)

source("inst/highdimregfunctions.R")

## convergence of marginal chain
load("output/highdimreg.meetings.RData")
niterations <- 2e3
ubounds <- sapply(1:niterations, function(t) tv_upper_bound(unlist(meetingtimes), lag, t))
g_tvbounds <- qplot(x = 1:niterations, y = ubounds, geom = "line") +
  ylab("TV upper bound") + xlab("iteration")
g_tvbounds <- g_tvbounds + geom_rangeframe() + ylim(0,1.1) + scale_y_continuous(breaks = c(0,1/2,1))
g_tvbounds
ggsave(filename = 'output/highdimreg.tvbounds.pdf', plot = g_tvbounds, width = 5, height = 4)

## trace and histogram
load(file = "output/highdimreg.longrun.RData")
gtrace <- ggplot(history.df %>% filter(iteration >= 10, iteration <= 1000, chain <= 3),
       aes(x = iteration, y = value, group = chain, colour = factor(chain))) + geom_line() + xlab("iteration") + ylab(TeX("$\\beta_{2564}$")) +
  scale_color_manual(values = c("black", "#005BBB", "#FFD500", "red")) + theme(legend.position = "none")
gtrace
ggsave(filename = 'output/highdimreg.trace.pdf', plot = gtrace, width = 5, height = 4)

burnin <- 1000
ghist <- ggplot(history.df %>% filter(iteration >= burnin), aes(x = value)) + geom_histogram(aes(y=..density..)) + xlab(TeX("$\\beta_{2564}$"))
ghist
ggsave(filename = 'output/highdimreg.hist.pdf', plot = ghist, width = 5, height = 4)

## unbiased asymptotic variance estimation
load("output/highdimreg.uavar.RData")


results.df <- foreach(irep = 1:nrep, .combine=rbind) %do% {
  run <- results[[irep]]
  data.frame(natoms = natoms_seq,
             rep = irep,
             estimator = run$estimator,
             cost = run$cost,
             pih = mean(run$pih),
             varh = run$varh,
             cost_fishyestimation = run$cost_fishyterms,
             fishyterms = run$fishyterms)
}
head(results.df)  
tail(results.df)  


max(results.df$cost)
max(results.df$cost-results.df$cost_fishyestimation)/2

table <- results.df %>% group_by(natoms) %>% summarise(
  estimate = mean(estimator),
  twostderror = 2 * sqrt(var(estimator)/nrep),
  totalcost = mean(cost),
  fishycost = mean(cost_fishyestimation),
  variance = var(estimator),
  variancepretty = prettyNum(var(estimator), digits = 2, scientific=T))
table$inefficiency <- prettyNum(table$totalcost * table$variance, digits=2, scientific =T)
table <- table %>% select(-variance) %>% setNames(c("R", "estimate  $v^{\\parallel}(P,\\test)$", "$2\\times \\hat{\\sigma}$", "cost", "fishy cost", 
                                                    "variance of $\\hat{v}(P,\\test)$", "inefficiency"))
print(table)
knitr::kable(table, digits = 1, row.names = NA, format = 'latex', escape = FALSE) %>%
  cat(., file = 'output/highdimreg.summary.tex')

cat("average run (in seconds):", mean(sapply(results, function(x) x$elapsedtime)), "\n")
cat("longest run (in seconds):", max(sapply(results, function(x) x$elapsedtime)), "\n")
hist(sapply(results, function(x) x$cost))
## compare unbiased MCMC efficiency with long-run MCMC efficiency
umcmc.efficiency.df <- foreach(irep = 1:nrep, .combine=rbind) %do% {
  cost = results[[irep]]$cost[1] - results[[irep]]$cost_fishyterms[1]
  pih = results[[irep]]$pih
  data.frame(rep = irep, pih = pih, cost = cost)
}
head(umcmc.efficiency.df)
umcmc.efficiency.summary.df <- umcmc.efficiency.df %>% summarise(v = 2*var(pih), meancost = mean(cost)/2) %>% mutate(inefficiency = v * meancost)

umcmc.efficiency.summary.df$mcmcinefficiency <- tail(table$`estimate  $v^{\\parallel}(P,\\test)$`, 1)
umcmc.efficiency.summary.df %>% mutate(relativeinefficiency =  inefficiency / mcmcinefficiency)

