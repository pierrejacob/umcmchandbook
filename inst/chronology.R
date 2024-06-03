library(umcmchandbook)
library(tidyverse)
setmytheme()
library(latex2exp)
set.seed(1)
library(doParallel)
library(doRNG)
library(dplyr)
# register parallel cores
registerDoParallel(cores = detectCores()-2)

durations <- rgamma(100000, shape = 20, rate = 10)
hist(durations, xlim = c(0,5))

create_chrono_onerun <- function(totalbudget){
  df <- data.frame()
  durs <- c()
  while (totalbudget > 0){
    duration <- rgamma(1, shape = 20, rate = 10)
    totalbudget <- totalbudget - duration
    durs <- c(durs, duration)
  }
  starttime <- c(0, cumsum(durs))[1:(length(durs))]
  endtime <- cumsum(durs)
  return(data.frame(starttime = starttime, endtime = endtime, isample = 1:(length(durs))))
}

nrep <- 8

chronodf <- foreach(irep = 1:nrep, .combine = rbind) %dopar% {
  df <- create_chrono_onerun(10)
  df$rep <- irep
  df
}



g <- ggplot(chronodf, aes(y = rep, yend = rep, x = starttime+0.01, xend = endtime-0.01, colour = factor(isample %% 2))) + 
  geom_segment(lineend = "round", linejoin = "mitre", size = 4) + geom_point()
g <- g + theme(legend.position = "none") + scale_color_manual(values = c("black", "grey"))
g <- g + scale_y_continuous(breaks = 1:nrep) + ylab("Machines")
g <- g + scale_x_continuous(breaks = c()) + xlab("Time")
g


# ggsave(filename = "chronology.pdf", plot = g, width = 8, height = 6)


g <- ggplot(chronodf, aes(y = rep, yend = rep, x = starttime+0.01, xend = endtime-0.01)) + 
  geom_segment(lineend = "round", linejoin = "mitre", size = 1) + geom_point(size = 5)
g <- g + theme(legend.position = "none") + scale_color_manual(values = c("black", "grey"))
g <- g + scale_y_continuous(breaks = 1:nrep) + ylab("Machines")
g <- g + scale_x_continuous(breaks = c()) + xlab("Time")
g

ggsave(filename = "chronology.pdf", plot = g, width = 8, height = 6)



