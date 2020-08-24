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
library(ggbeeswarm)
library(ltmle)
library(Plasmode)
library(SuperLearner)
library(ctmle)
library(sl3)
library(tmle3)
library(tictoc)
library(doParallel)
library(doRNG)
library(foreach)
library(randomForest)

# require(purrr)
# require(furrr)
setwd("~/Desktop/HuangGroup/cvtmle_plasmode/Code")
set.seed(42782)
options(tibble.print_max = 40, tibble.print_min = 30)
registerDoParallel(cores=detectCores(all.tests = T)-2)

# out_path <- "D:/SICS/Projects/TMLE_Plasmode/"
# data_path <- "D:/SICS/Data and Instruments/Archive/"


# Set simulation parameters
{
  sims.ver <- "plas"
  # sims.ver <- "5var"
  # sims.ver <- "5var.then.plas"
  
  Effect_Size <- 6.6
  
  plas.seed <- 2222
  
  ########
  # parameters for plasmode
  #######
  data.ver <- "FULL"
  # data.ver <- "13"
  size <- 1178
  # size <-200
  plas_sim_N <- 500; use.subset <- F
  generateA <- F
  

  
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
}else{
  sim_boots <- sims.obj
}

# DC implementation
source("20200705-DCDR-Functions.R")
# getRES function is now relocated to a separate file (20200720-Algos-code.R)
source("20200720-Algos-code.R")

##################################
## PARALLELIZE ANALYSES
##################################

{
set.seed(42782)
tic()
N_sims <- 100 # this should <= plas_sim_N
# regression models for GLM / AIPW
reg.formulas <- make.formula("Y", "A",ver = data.ver,sims.ver = sims.ver)
expForm <- reg.formulas$expForm
outForm <- reg.formulas$outForm

# specify which set of learners for SL
# ### NON-SMOOTH
# short_tmle_lib <- SL_param_list
# tmle_lib <- lrnr_SL
# aipw_lib <- SL.lib
# 
# SMOOTH
short_tmle_lib <- SL_list
tmle_lib <- lrnr_SL_param
aipw_lib <- SL.param


# set .errorhandling="remove" if want to discard
# set .errorhandling="stop" by default

# boot1 <- foreach(i = 1:N_sims,.errorhandling="stop") %dopar% {
boot1 <- foreach(i = 1:N_sims,.errorhandling="remove") %dopar% {
  require(tidyverse)
  require(tmle3)
  require(sl3)
  require(SuperLearner)

  # Initialize dataset
  # i<-10 # DEBUG
  if (sims.ver == "plas" | sims.ver =="5var.then.plas"){
    # i <- 1
    plas_data <- cbind(id = plas_sims$Sim_Data[i],
                       A = plas_sims$Sim_Data[i + (2*plas_sim_N)],
                       Y = plas_sims$Sim_Data[i + plas_sim_N])
    colnames(plas_data) <- c("id", "A", "Y")
    set1 <- suppressMessages(left_join(as_tibble(plas_data), as_tibble(plas))) #dplyr::select(as_tibble(plas), -Y5, -A1))) # add covariates
    tset <- set1 %>% mutate(YT = (Y-min(set1$Y))/(max(set1$Y)- min(set1$Y)))
    
  }else{
    # Initialize dataset
    # i<-2
    ss <- 600
    set1 <- as_tibble(cbind(C1 = sim_boots[[i]]$C1[1:ss],
                            C2 = sim_boots[[i]]$C2[1:ss],
                            C3 = sim_boots[[i]]$C3[1:ss],
                            C4 = sim_boots[[i]]$C4[1:ss],
                            C5 = sim_boots[[i]]$C5[1:ss],
                            A = sim_boots[[i]]$A[1:ss],
                            Y = sim_boots[[i]]$Y[1:ss]))
    tset <- set1 %>% mutate(YT = (Y-min(set1$Y))/(max(set1$Y)- min(set1$Y))) # generate a bounded Y for TMLE
  }
  getRES(set1, tset, aipw_lib, tmle_lib, short_tmle_lib,
         doIPW = 1,
         doAIPW=0, doDCAIPW=0,
         doManuTMLE=0, doShortTMLE = 0,
         doDCTMLE=0,
         num_cf=5,
         #control=list()
         control=SuperLearner.CV.control(V=2)
         )
}
toc()
}

(N_boot <- length(boot1))



##################################
## SUMMARIZE AND VISUALIZE
##################################
{
source("20200816-Result-Summary.R")
cat("plas.seed=",plas.seed, "generateA=",generateA, "sims.ver=",sims.ver)
summarise.res(boot1)
}


##################################
## DEBUG
##################################
# {
# # append results
# meta_res_bk <- rbind(meta_res_bk, 
#                   cbind(N = nrow(plas), iter = N_boot, 
#                         A_mu_ATE = sim_res$mu_ATE[1], A_med_ATE = sim_res$med_ATE[1],
#                         A_mu_SE = sim_res$mu_SE[1], A_med_SE = sim_res$med_ATE[1],
#                         A_mu_bias = sim_res$mu_bias[1], A_med_bias = sim_res$med_bias[1],
#                         A_MSE = sim_res$MSE[1], A_coverage = sim_res$coverage[1],
#                         T_mu_ATE = sim_res$mu_ATE[2], T_med_ATE = sim_res$med_ATE[2],
#                         T_mu_SE = sim_res$mu_SE[2], T_med_SE = sim_res$med_ATE[2],
#                         T_mu_bias = sim_res$mu_bias[2], T_med_bias = sim_res$med_bias[2],
#                         T_MSE = sim_res$MSE[2], T_coverage = sim_res$coverage[2],
#                         aipw_lib = toString(aipw_lib), tmle_lib = tmle_lib$name))
# }

#writexl::write_xlsx(data.frame(meta_res), paste0(out_path,"sim_results_v1.xlsx"))

# writeClipboard(meta_res_bk[15,])
# writeClipboard(meta_res_bk[2,])

# example_plasmode <- set1
# writexl::write_xlsx(example_plasmode, paste0(out_path, "plasmode_data_200.xlsx"))
  
# meta_res_bk <- meta_res
# res1 <- c(meta_res[seq(1,length(meta_res),2)])
# res2 <- c(meta_res[seq(2,length(meta_res),2)])
# writeClipboard(as.character(cbind(res1)))
# writeClipboard(as.character(cbind(res2)))







##################################
##################################
# ## diagnostic for SL-estimated PS
# set1 <- set1 %>% add_column(tmle_ps = model_tmle$cum.g[,,1], tmle_ps_untrunc = model_tmle$cum.g.unbounded[,,1])
# ggplot(data = set1) + 
#   geom_histogram(aes(x = denom, 
#   #geom_histogram(aes(x = tmle_ps_untrunc, 
#                      fill = as.factor(A1)), alpha = 0.4, position = "dodge", binwidth = 0.05) +
#   labs(title = "Propensity scores estimated by SL (TMLE), by treatment status",
#        fill = "A1" ) + scale_x_continuous(breaks = seq(0,10,0.1)) +
#   theme(legend.position = c(0.1,0.9), legend.background = element_blank())
# 
# set1 %>% group_by(A1) %>% summarize(N = n(), min(tmle_ps), mean(tmle_ps), max(tmle_ps))
# set1 %>% group_by(A1) %>% summarize(N = n(), min(denom), mean(denom), max(denom))
# set1 %>% group_by(A1) %>% summarize(N = n(), min(wt), mean(wt), max(wt))
# 
# # use % truncated as a measure of overfit (?)
# set1 %>% group_by(A1) %>% 
#   summarize(sum(tmle_ps == 0.1), mean(tmle_ps == 0.1),
#             sum(tmle_ps == 0.9), mean(tmle_ps == 0.9)) 
# 
# set1 %>% dplyr::select(A1, tmle_ps_untrunc, denom) %>% gather(model, val, -A1) %>% 
# ggplot() + geom_density(aes(val, fill = factor(A1)), alpha = 0.4) +
#   labs(title = "Propensity scores estimated by SL (TMLE, truncated at 0.9), by treatment status",
#        fill = "A1" ) + scale_x_continuous(breaks = seq(0,10,0.1)) +
#   theme(legend.position = c(0.1,0.9), legend.background = element_blank()) +
#   facet_wrap(~model)
