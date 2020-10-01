####################################################################
# Compare performance of CV-TMLE vs. AIPW vs. IPW-GLM in realistic sample size
# PLASMODE SIMULATION
#
# Created: 2020 9 March
# Modified: 2020 9 June
#
# Updated with fixed plasmode simulation function (PlasmodeContNew)
#
####################################################################

library(tidyverse)
library(ggplot2)
library(haven)
library(readxl)
library(table1)
# library(ggbeeswarm)
# library(ltmle)
library(Plasmode)
library(SuperLearner)
#library(ctmle)
library(sl3)
library(tmle3)
library(tictoc)
library(doParallel)
library(doRNG)
library(foreach)
library(randomForest)


# path <- "~/Desktop/HuangGroup/cvtmle_plasmode"
# setwd(paste0(path,"/Code"))
path <- "/n/holyscratch01/murphy_lab/Users/xmeng/submitted092920"
setwd(paste0(path))
set.seed(42782)
options(tibble.print_max = 40, tibble.print_min = 30)
no_cores <- detectCores(all.tests = T) - 2
registerDoParallel(cores=no_cores)


# Set simulation parameters
{
  sims.ver <- "plas"
  # sims.ver <- "5var"
  # sims.ver <- "5var.then.plas"
  # sims.ver <- "5var.hard"
  
  Effect_Size <- 6.6
  
  plas.seed <- 1111
  
  ########
  # parameters for plasmode
  #######
  # p=331 for the whole set 
  data.ver <- "FULL"
  # data.ver <- "13"
  size <- 1178
  # size <-200
  plas_sim_N <- 500; use.subset <- F
  generateA <- T
  estimateWithMore <- F; p.sim = 100;p.est <- 50
  randVar = F # Do we permute variable order?
  isHandPick = T; idx.handpick <- c(1,2,5,18, 217)
  interact=T
  
  
  ########
  # parameters for 5 var, 5var.then.plas
  #######
  # Nsets <- 10000
  # Nsamp <- 3000
  Nsets <- 500
  Nsamp <- 600
}

source("20200803-Sims-Function.R")
cat(sims.ver)
sims.obj <- general.sim(sims.ver)

if (sims.ver == "plas"|sims.ver == "5var.then.plas"){
  plas <- sims.obj$plas
  plas_sims <- sims.obj$plas_sims
  vars <- sims.obj$vars
}else{
  sim_boots <- sims.obj
}

# Plot the regression coefficient in Plasmode simulation
# plot(plas_sims$TrueOutBeta)


##################################
## PARALLELIZE ANALYSES
##################################
# Set param
{
  # DC implementation
  source("20200705-DCDR-Functions.R")
  # getRES function is now relocated to a separate file (20200720-Algos-code.R)
  source("20200720-Algos-code.R")


non.par <- T
# specify which set of learners for SL
if (non.par == T){
  ### NON-SMOOTH
  short_tmle_lib <- SL_param_list
  tmle_lib <- lrnr_SL
  aipw_lib <- SL.lib
}
else{
  # SMOOTH
  short_tmle_lib <- SL_list
  tmle_lib <- lrnr_SL_param
  aipw_lib <- SL.param
}


# aipw_lib <- c("SL.glmnet")

# errorhandling="stop"
errorhandling="remove"

est.interact <- F # Estimate with first order interaction?

N_sims <- 50# this should <= plas_sim_N

doIPW = 0; doLASSO=0;
doAIPW=0; doDCAIPW=1
doManuTMLE=0; doDCTMLE=0
num_cf=5
doShortTMLE = 0
#control=list()
control=SuperLearner.CV.control(V=2)
}


#########
# run
##########
source("20200904-run-sim-code.R")

# print(paste("It takes",round((ptm1 - ptm0)[3],2),"s"))
print(paste("Cores used:",no_cores))
print(paste("plas.seed=",plas.seed, "generateA=",generateA, "sims.ver=",sims.ver))

save(boot1, file=paste("./RDataFiles/result.RData",sep=""))  





