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
library(foreach)
library(randomForest)

set.seed(42782)
options(tibble.print_max = 40, tibble.print_min = 30)
registerDoParallel(cores = detectCores(all.tests = T))

# out_path <- "D:/SICS/Projects/TMLE_Plasmode/"
# data_path <- "D:/SICS/Data and Instruments/Archive/"

##########################################
# Plasmode simulation
##########################################
# org_data <- read_dta(paste0(data_path,"merged_MI.dta")) %>% filter(visit == 0) 
# vars <- c(names(org_data)[12:21], names(org_data)[26:41], names(org_data)[61:83], names(org_data)[140:421])
# org_data_imp <- org_data %>% mutate(A1 = as.double(ppBMI > 27)) %>% select(Y5 = bmi, A1, vars) %>% 
#   mutate_if(is.character, as.factor) %>% mutate_all(list(~as.double(.))) %>%
#   mutate_all(list(~if_else(is.na(.), mean(., na.rm = T), .)))  
# plas <- org_data_imp %>% add_column(id = c(1:nrow(org_data_imp))) %>% mutate(id = paste0("ID",id)) %>% data.frame(.)
# varnames <- paste0(rep("VAR_", length(vars)), as.character(seq(1,length(vars)), 1))
# plas %>% rename_at(vars(vars), ~varnames) %>% write_dta(., paste0(out_path,"plas_data.dta"))

plas <- haven::read_dta(paste0(out_path,"plas_data.dta")) 
vars <- names(plas[3:333])

### Outcome and exposure with high-dim covars (N = 331 covars)
  outForm <- "Y5 ~ A1"
  expForm <- "A1 ~ "
  for(i in 1:length(vars)){
    outForm <- paste0(outForm," + ", vars[i])
    if(i == 1){expForm <- paste0(expForm, vars[i])}
    else{expForm <- paste0(expForm, " + ", vars[i])}
  }
 
### Outcome and exposure model with original, restricted covariate set (N = 16 covars)
  # plas_org <- plas %>%
  #   dplyr::select(A1,
  #                 L0.a = VAR_12, L0.b = VAR_4, L0.c = VAR_2, L0.d = VAR_11, L0.e = VAR_217, L0.f = VAR_13, L0.g = VAR_1,
  #                 L0.h = VAR_49, L0.i = VAR_42, L0.j = VAR_41, L0.k = VAR_3,
  #                 L1.a = VAR_9, L1.b = VAR_27, L1.c = VAR_28, L1.d = VAR_29,
  #                 L1.e = VAR_31, L1.f = VAR_16, L1.g = VAR_8, Y5)
  # plas_org <- plas_org %>% add_column(id = c(1:nrow(plas_org))) %>% mutate(id = paste0("ID",id)) %>% data.frame(.)
  # outForm <- Y5 ~ A1 + L0.a + L0.b + L0.c + L0.d + L0.e + L0.f + L0.g + L0.h + L0.i + L0.j + L0.k
  # expForm <- A1 ~ L0.a + L0.b + L0.c + L0.d + L0.e + L0.f + L0.g + L0.h + L0.i + L0.j + L0.k
  
plas <- data.frame(plas)

plas_sim_N <- 500
Effect_Size <- 6.6 # simulated risk difference = large change (e.g. absolute units)
#Effect_Size <- 0.66 # simulated risk difference = 2/3 SD higher BMI
plas_sims <- PlasmodeContNew(formulaOut = as.formula(outForm), objectOut = NULL,
                             formulaExp = as.formula(expForm), objectExp = NULL,
                             data = plas, idVar = "id", 
                             effectOR = Effect_Size, MMOut = 1, MMExp = 1, 
                             nsim = plas_sim_N, size = nrow(plas),
                             exposedPrev = NULL)

##########################################
# Initialize necessary parameters for estimator
##########################################
# TMLE parameters
#SL.lib <- c("SL.randomForest", "SL.xgboost", "SL.nnet", "SL.glm", "SL.glmnet", "SL.polymars")
#SL.lib <- list("SL.randomForest", "SL.xgboost", "SL.nnet", "SL.glm", c("SL.glmnet", "All"), c("SL.polymars", "All"))
SL.lib <- list("SL.randomForest", "SL.xgboost", "SL.glm", c("SL.polymars", "All"))
SL.lib.tmle <- c("SL.randomForest", "SL.xgboost", "SL.glm", "SL.polymars")
SL.param <- c("SL.glm", "SL.glmnet", "SL.polymars")

# Specify the NP-SEM for the TMLE - including bounded, tranformed Y ("YT")
npsem <- list(define_node("Z", vars),
                          #c("L0.a", "L0.b", "L0.c", "L0.d", "L0.e", "L0.f", "L0.g", "L0.h", "L0.i", "L0.j", "L0.k")),
              define_node("A", c("A"), c("Z")),
              define_node("Y", c("YT"), c("A", "Z")))

# Specify the learners for CV-TMLE
lrnr_SL <- make_learner(Lrnr_pkg_SuperLearner, SL.lib.tmle)
lrnr_SL_param <- make_learner(Lrnr_pkg_SuperLearner, SL.param)
lrnr_glm <- make_learner(Lrnr_glm_fast)
lrnr_mean <- make_learner(Lrnr_mean)

SL_list <- list(Y = lrnr_SL, A = lrnr_SL)
SL_param_list <- list(Y = lrnr_SL_param, A = lrnr_SL_param)
glm_list <- list(Y = lrnr_glm_fast, A = lrnr_glm_fast)

# regression models for GLM / AIPW
#expForm <- "A ~ L0.a + L0.b + L0.c + L0.d + L0.e + L0.f + L0.g + L0.h + L0.i + L0.j + L0.k"
#outForm <- "Y ~ A + L0.a + L0.b + L0.c + L0.d + L0.e + L0.f + L0.g + L0.h + L0.i + L0.j + L0.k"

###############################
# FUNCTION TO RUN AIPW, CV-TMLE, and IPW+GLM
###############################
getRES <- function(gset, tset, aipw_lib = NULL, tmle_lib = NULL, doAIPW = 0, doIPW = 0, doTMLE = 0){
  res <- NULL
  
  # AIPW
  if (doAIPW == 1){
  outcome <- as.character(as.formula(outForm)[[2]])
  exposure <- as.character(as.formula(expForm)[[2]]) 
  out_vec <- pull(select(gset, outcome))
  exp_vec <- pull(select(gset, exposure))
  
  denom_fit <- glm(data = gset, formula = as.formula(expForm), family = "binomial")
  gmodel <- glm(data = gset, formula = as.formula(outForm), family = "gaussian")
  
  expVars <- names(denom_fit$coefficients)[2:length(denom_fit$coefficients)]
  exp_pred_mat <- select(gset, expVars)
  suppressMessages(learnPS <- SuperLearner(Y = exp_vec, X = exp_pred_mat, SL.library = aipw_lib, family = binomial()))
  pred_ps <- c(learnPS$SL.predict)
    
  outVars <- names(gmodel$coefficients)[2:length(gmodel$coefficients)]
  out_pred_mat <- select(gset, outVars)
  out_x1_mat <- out_x0_mat <- out_pred_mat
  out_x1_mat[exposure] <- rep(1, length(exp_vec))
  out_x0_mat[exposure] <- rep(0, length(exp_vec))
    
  learnOut_1 <- SuperLearner(Y = out_vec, X = out_pred_mat, newX = out_x1_mat, SL.library = aipw_lib)
  pred_x1 <- c(learnOut_1$SL.predict)
  learnOut_0 <- SuperLearner(Y = out_vec, X = out_pred_mat, newX = out_x0_mat, SL.library = aipw_lib)
  pred_x0 <- c(learnOut_0$SL.predict)
  
  gset <- gset %>% mutate(ps_u = 1-pred_ps, ps_t = pred_ps, pred_u = pred_x0, pred_t = pred_x1)
  # ps_u = (1-denom_fit$fitted.values), ps_t = denom_fit$fitted.values,
  # pred_u = predict(gmodel, data.frame(C1 = gset$C1, C2 = gset$C2, C3 = gset$C3,
  #                                     C4 = gset$C4, C5 = gset$C5, A = rep(0, nrow(gset)))),
  # pred_t = predict(gmodel, data.frame(C1 = gset$C1, C2 = gset$C2, C3 = gset$C3,
  #                                     C4 = gset$C4, C5 = gset$C5, A = rep(1, nrow(gset)))))
  aipw <- RCAL::ate.aipw(y = gset$Y, tr = gset$A, mfp = cbind(gset$ps_u, gset$ps_t), mfo = cbind(gset$pred_u, gset$pred_t))
  ATE <- aipw$diff.est[2]
  SE <- sqrt(aipw$diff.var[2])
  TYPE <- "AIPW"
  res <- rbind(res, cbind(ATE, SE, TYPE))
  }
  
  # IPW + GLM
  if (doIPW == 1){
  num <- summary(glm(data = gset, formula = A ~ 1, family = "gaussian"))$coefficient[1,1]
  denom_fit <- glm(data = gset,
                   formula = as.formula(expForm),
                   family = "binomial")
  gset <- gset %>% add_column(denom = denom_fit$fitted.values) %>%
    mutate(wt = if_else(A == 1, num/denom, (1-num)/(1-denom))) #Stabilized IPTW
  gset <- gset %>% mutate(wt2 = case_when(wt > quantile(gset$wt, 0.95) ~ quantile(gset$wt, 0.95),
                                          wt < quantile(gset$wt, 0.05) ~ quantile(gset$wt, 0.05),
                                          T ~ wt)) # Truncated IPTW
  model_glm <- glm(data = gset, weight = wt2,
                   formula = as.formula(outForm),
                   family = "gaussian")
  ATE <- summary(model_glm)$coefficients[2,1]
  SE <- summary(model_glm)$coefficients[2,2]
  TYPE <- "GLM + IPW"
  res <- rbind(res, cbind(ATE, SE, TYPE))
  }
  
  # CV-TMLE default specification
  # nodes <- list(W = c("C1","C2","C3","C4","C5"), A = "A", Y = "Y")
  # tmle3_autofit <- tmle3(tmle_ATE(treatment_level = 1, control_level = 0), data.table::copy(tset), nodes, learner_list)
  # ATE <- tmle3_autofit$summary$tmle_est
  # SE <- tmle3_autofit$summary$se
  
  # CV-TMLE full specification
  if (doTMLE == 1){
  tmle_task <- tmle3_Task$new(tset, npsem = npsem)
  factor_list <- list(define_lf(LF_emp, "Z"), define_lf(LF_fit, "A", learner = tmle_lib), define_lf(LF_fit, "Y", learner = tmle_lib))
  likelihood_def <- Likelihood$new(factor_list)
  likelihood <- likelihood_def$train(tmle_task)
  ate_params <- list(Param_ATE$new(likelihood, define_lf(LF_static, "A", value = 1), define_lf(LF_static, "A", value = 0)))
  updater <- tmle3_Update$new()
  targeted_likelihood <- Targeted_Likelihood$new(likelihood, updater)
  tmle3_fit <- fit_tmle3(tmle_task, targeted_likelihood, ate_params, updater)
  ATE <- tmle3_fit$summary$tmle_est * (max(tset$Y)-min(tset$Y)) # backtransform ATE
  SE <- tmle3_fit$summary$se * (max(tset$Y)-min(tset$Y)) # backtransform SE
  TYPE <- "CV-TMLE"
  res <- rbind(res, cbind(ATE, SE, TYPE))
  }
  
  return(res)
}


#meta_res <- NULL # intialize the results table

##################################
## PARALLELIZE ANALYSES
##################################
{
set.seed(42782)
tic()
N_sims <- 100

# specify which set of learners for SL
#### NON-SMOOTH
# tmle_lib <- lrnr_SL
# aipw_lib <- SL.lib

#### SMOOTH
tmle_lib <- lrnr_SL_param
aipw_lib <- SL.param


boot1 <- foreach(i = 1:N_sims) %dopar% {
  require(tidyverse)
  require(tmle3)
  require(sl3)
  require(SuperLearner)

  # Initialize dataset
  plas_data <- cbind(id = plas_sims$Sim_Data[i],
                     A = plas_sims$Sim_Data[i + (2*plas_sim_N)],
                     Y = plas_sims$Sim_Data[i + plas_sim_N])
  colnames(plas_data) <- c("id", "A", "Y")
  set1 <- suppressMessages(left_join(as_tibble(plas_data), as_tibble(plas))) #dplyr::select(as_tibble(plas), -Y5, -A1))) # add covariates
  tset <- set1 %>% mutate(YT = (Y-min(set1$Y))/(max(set1$Y)- min(set1$Y)))
  
  getRES(set1, tset, aipw_lib, tmle_lib, doAIPW = 1)
}
toc()
}

#meta_res_bk <- meta_res
#meta_res_bk <- NULL

##################################
## SUMMARIZE AND VISUALIZE
##################################
{
tmp <- NULL
for(i in 1:N_sims){
  tmp <- rbind(tmp, boot1[[i]])
  }

sim_corr1 <- as_tibble(tmp) %>% 
  mutate(ATE = as.double(ATE), SE = as.double(SE),
         lb = as.double(ATE) - 1.96*as.double(SE), 
         ub = as.double(ATE) + 1.96*as.double(SE)) %>% 
  add_column(iter = c(1:nrow(tmp))) %>% mutate(bias = as.double(ATE) - Effect_Size)

#write_csv(plas_corr1, paste0(out_path, Sys.Date(), "-plasmode_sim_CV-TMLE_GLM_IPW_1000.csv"))

sim_res <- sim_corr1 %>% group_by(TYPE) %>% 
  summarize(mu_ATE = mean(ATE), med_ATE = median(ATE), 
            mu_SE = mean(SE), med_SE = median(SE), 
            mu_bias = mean(bias), med_bias = median(bias), 
            var = var(ATE), MSE = var + mu_bias^2,
            coverage = sum(lb <= Effect_Size & ub >= Effect_Size)/N_sims)

pos <- max(sim_corr1$ub)
sim_plot <- sim_corr1 %>% 
  ggplot(aes(x = iter/2, y = ATE, color = TYPE)) +
  geom_point() + geom_errorbar(aes(ymin = lb, ymax = ub)) + 
  geom_hline(aes(yintercept = Effect_Size)) + 
  # geom_hline(data = filter(sim_corr1, substring(TYPE,1,3) %in% c("AIP", "GLM")), aes(yintercept = mean(ATE)), linetype = "dashed") +
  geom_hline(data = filter(sim_corr1, substring(TYPE,1,2) == "CV"), aes(yintercept = mean(ATE)), linetype = "dotted") +
  # geom_text(data = filter(sim_corr1, substring(TYPE,1,3) %in% c("AIP", "GLM")),
  #           aes(N_sims, pos, hjust = "right",
  #               label = paste("\n\nMean bias = ", as.character(round(mean(bias), 3)),
  #                             "\nMedian bias = ", as.character(round(median(bias), 3)),
  #                             "\nMSE = ", as.character(round(var(ATE) + mean(bias)^2, 3)),
  #                             "\nCoverage = ", as.character(round(sum(lb <= Effect_Size & ub >= Effect_Size)/N_sims, 3)) ))) +
  geom_text(data = filter(sim_corr1, substring(TYPE,1,2) == "CV"),
            aes(N_sims, pos, hjust = "right",
                label = paste("\n\nMean bias = ", as.character(round(mean(bias), 3)),
                              "\nMedian bias = ", as.character(round(median(bias), 3)),
                              "\nMSE = ", as.character(round(var(ATE) + mean(bias)^2, 3)),
                              "\nCoverage = ", as.character(round(sum(lb <= Effect_Size & ub >= Effect_Size)/N_sims, 3)) ))) +
  labs(x = "iteration", color = "Estimator") +
  facet_wrap(~TYPE) +
  theme(legend.position = "none")
sim_plot

# append results
meta_res_bk <- rbind(meta_res_bk, 
                  cbind(N = nrow(plas), iter = N_sims, 
                        A_mu_ATE = sim_res$mu_ATE[1], A_med_ATE = sim_res$med_ATE[1],
                        A_mu_SE = sim_res$mu_SE[1], A_med_SE = sim_res$med_ATE[1],
                        A_mu_bias = sim_res$mu_bias[1], A_med_bias = sim_res$med_bias[1],
                        A_MSE = sim_res$MSE[1], A_coverage = sim_res$coverage[1],
                        T_mu_ATE = sim_res$mu_ATE[2], T_med_ATE = sim_res$med_ATE[2],
                        T_mu_SE = sim_res$mu_SE[2], T_med_SE = sim_res$med_ATE[2],
                        T_mu_bias = sim_res$mu_bias[2], T_med_bias = sim_res$med_bias[2],
                        T_MSE = sim_res$MSE[2], T_coverage = sim_res$coverage[2],
                        aipw_lib = toString(aipw_lib), tmle_lib = tmle_lib$name))
}

#writexl::write_xlsx(data.frame(meta_res), paste0(out_path,"sim_results_v1.xlsx"))

writeClipboard(meta_res_bk[15,])
writeClipboard(meta_res_bk[2,])

# example_plasmode <- set1
# writexl::write_xlsx(example_plasmode, paste0(out_path, "plasmode_data_200.xlsx"))
  
# meta_res_bk <- meta_res
# res1 <- c(meta_res[seq(1,length(meta_res),2)])
# res2 <- c(meta_res[seq(2,length(meta_res),2)])
# writeClipboard(as.character(cbind(res1)))
# writeClipboard(as.character(cbind(res2)))









##################################
##################################
## diagnostic for SL-estimated PS
set1 <- set1 %>% add_column(tmle_ps = model_tmle$cum.g[,,1], tmle_ps_untrunc = model_tmle$cum.g.unbounded[,,1])
ggplot(data = set1) + 
  geom_histogram(aes(x = denom, 
  #geom_histogram(aes(x = tmle_ps_untrunc, 
                     fill = as.factor(A1)), alpha = 0.4, position = "dodge", binwidth = 0.05) +
  labs(title = "Propensity scores estimated by SL (TMLE), by treatment status",
       fill = "A1" ) + scale_x_continuous(breaks = seq(0,10,0.1)) +
  theme(legend.position = c(0.1,0.9), legend.background = element_blank())

set1 %>% group_by(A1) %>% summarize(N = n(), min(tmle_ps), mean(tmle_ps), max(tmle_ps))
set1 %>% group_by(A1) %>% summarize(N = n(), min(denom), mean(denom), max(denom))
set1 %>% group_by(A1) %>% summarize(N = n(), min(wt), mean(wt), max(wt))

# use % truncated as a measure of overfit (?)
set1 %>% group_by(A1) %>% 
  summarize(sum(tmle_ps == 0.1), mean(tmle_ps == 0.1),
            sum(tmle_ps == 0.9), mean(tmle_ps == 0.9)) 

set1 %>% dplyr::select(A1, tmle_ps_untrunc, denom) %>% gather(model, val, -A1) %>% 
ggplot() + geom_density(aes(val, fill = factor(A1)), alpha = 0.4) +
  labs(title = "Propensity scores estimated by SL (TMLE, truncated at 0.9), by treatment status",
       fill = "A1" ) + scale_x_continuous(breaks = seq(0,10,0.1)) +
  theme(legend.position = c(0.1,0.9), legend.background = element_blank()) +
  facet_wrap(~model)
