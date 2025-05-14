#####################################################
# SET THE VALUES OF THE PARAMETERS THAT DONT CHANGE #
#####################################################
set.seed(8)

source("sim_datasets.R") # here the simulation of the data sets (Simulation_cristina) 
source("Simulation_function.R")  # here the simulation functions 

#####################################################
# SET THE VALUES OF THE PARAMETERS THAT DONT CHANGE #
#####################################################

#This is always fixed
n <- 500 # number of subjects
K <- 4 # number of planned repeated measurements per subject, per outcome
t.max <- 6.5 # maximum follow-up time
nu <- 1.75 # shape for the Weibull baseline hazard
gamma <- c(-1.6) 

cens <- 0.1 # the random censoring

#Fixed in scenario 1

betas <- c("Intercept" = 2.17 , "time" = -0.07) 
ntrial <- 8
D <- matrix(c(2.45,0,0,0.07),ncol=2)


#############
# SIMULATE  #
#############
nSims <- 500
out <- simulation.function(nSims = nSims,n = n,K=K,t.max=t.max,betas=betas,phi=phi,ntrial=ntrial,alpha=alpha,nu=nu,gamma=gamma,D=D,cens=cens) # get the results of the simulations for the given values

