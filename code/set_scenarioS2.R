rm(list=ls())

library("MASS")
library("PROreg")
library("splines")
library("statmod")
library("rstan")
set.seed(8)



phi<- 0.05
alpha<- -1.5

source("runSimsS2_lp.R")
save.image("S2_p005a_15.RData")