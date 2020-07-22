####################################################################
# Compare performance of CV-TMLE vs. AIPW vs. IPW-GLM in realistic sample size
# SIMPLE PARAMETRIC SIMULATION
#
# Created: 2020 4 March
# Modified: 2020 9 June
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

# require(purrr)
# require(furrr)

setwd("~/Desktop/HuangGroup/cvtmle_plasmode/Code")
set.seed(42782)
options(tibble.print_max = 40, tibble.print_min = 30)
registerDoParallel(cores = detectCores(all.tests = T) - 1)

##################################
# Generate basic sim dataset
##################################
set.seed(42782)

sim_boots <- NULL
# Nsets <- 10000
# Nsamp <- 3000
Nsets <- 1000
Nsamp <- 600
Effect_Size <- 6.6

draw_sims <- function(i){
  require(tidyverse)
  
  C1 <- rnorm(Nsamp, 0, 1)
  C2 <- rnorm(Nsamp, C1 + 2, 2)
  C3 <- rnorm(Nsamp, 2, abs(2*C2))
  C4 <- rnorm(Nsamp, C2^2 + 2*C3, abs(C1))
  C5 <- rnorm(Nsamp, C3*C4, abs(C2-C1))
  as_tibble(cbind(C1, C2, C3, C4, C5)) %>% summarize_all(list(~mean(.)))
  Pr_A = plogis(C1 + C2/20 + C3/50 + C4/200 + C5/5000)
  A <- rbinom(length(Pr_A), 1, Pr_A)
  Y <- Effect_Size*A + 10*C1 + 0.5*C2^2 + 0.66*C3 + 0.25*C4 + 0.01*C3*C4 + 4*log(C5^2)
  I <- i
  res <- c(as_tibble(cbind(C1, C2, C3, C4, C5, A, Y, I)))
  return(res)
}

sim_boots <- foreach(i = 1:Nsets) %dopar% {
  draw_sims(i)
}



sim_boots2 <- sim_boots[1:1000]



##########################################
# Initialize necessary parameters for estimator
##########################################
{# TMLE parameters
#SL.lib <- c("SL.randomForest", "SL.xgboost", "SL.nnet", "SL.glm", "SL.glmnet", "SL.polymars")
#SL.lib <- list("SL.randomForest", "SL.xgboost", "SL.nnet", "SL.glm", c("SL.glmnet", "All"), c("SL.polymars", "All"))
SL.lib <- list("SL.randomForest", "SL.xgboost", "SL.glm", c("SL.polymars", "All"))
SL.lib.tmle <- c("SL.randomForest", "SL.xgboost", "SL.glm", "SL.polymars")
SL.param <- c("SL.glm", "SL.glmnet", "SL.polymars")

# Specify the NP-SEM for the TMLE - including bounded, tranformed Y ("YT")
npsem <- list(define_node("Z", 
                          c("C1", "C2", "C3", "C4", "C5")),
              define_node("A", c("A"), c("Z")),
              define_node("Y", c("YT"), c("A", "Z")))

# Specify the learners for CV-TMLE
lrnr_SL <- make_learner(Lrnr_pkg_SuperLearner, SL.lib.tmle)
lrnr_SL_param <- make_learner(Lrnr_pkg_SuperLearner, SL.param)
lrnr_glm <- make_learner(Lrnr_glm_fast)
lrnr_mean <- make_learner(Lrnr_mean)

SL_list <- list(Y = lrnr_SL, A = lrnr_SL)
SL_param_list <- list(Y = lrnr_SL_param, A = lrnr_SL_param)
glm_list <- list(Y = lrnr_glm, A = lrnr_glm)
}


# ###############################
# # FUNCTION TO RUN AIPW, CV-TMLE, and IPW+GLM
# ###############################
# getRES <- function(gset, tset, aipw_lib = NULL, tmle_lib = NULL){
#   # # AIPW
#   # outcome <- as.character(as.formula(outForm)[[2]])
#   # exposure <- as.character(as.formula(expForm)[[2]])
#   # out_vec <- pull(select(gset, outcome))
#   # exp_vec <- pull(select(gset, exposure))
#   # 
#   # # denom_fit <- glm(data = gset, formula = as.formula(expForm), family = "binomial",fit=FALSE)
#   # # gmodel <- glm(data = gset, formula = as.formula(outForm), family = "gaussian",fit=FALSE)
#   # #
#   # # expVars <- names(denom_fit$coefficients)[2:length(denom_fit$coefficients)]
#   # # exp_pred_mat <- data.frame(select(gset, expVars))
#   # exp_pred_mat <- model.matrix(object = as.formula(expForm), data=gset)[,-1]
#   # exp_pred_mat <- data.frame(exp_pred_mat)
#   # suppressMessages(learnPS <- SuperLearner(Y = exp_vec, X = exp_pred_mat, SL.library = aipw_lib, family = binomial()))
#   # pred_ps <- c(learnPS$SL.predict)
#   # 
#   # # outVars <- names(gmodel$coefficients)[2:length(gmodel$coefficients)]
#   # # out_pred_mat <- data.frame(select(gset, outVars))
#   # out_pred_mat <- model.matrix(object = as.formula(outForm), data=gset)[,-1]
#   # out_pred_mat <- data.frame(out_pred_mat)
#   # out_x1_mat <- out_x0_mat <- out_pred_mat
#   # out_x1_mat[exposure] <- rep(1, length(exp_vec))
#   # out_x0_mat[exposure] <- rep(0, length(exp_vec))
#   # 
#   # # learnOut_1 <- SuperLearner(Y = as.matrix(out_vec), X = out_pred_mat, newX = out_x1_mat, SL.library = aipw_lib)
#   # # pred_x1 <- c(learnOut_1$SL.predict)
#   # # learnOut_0 <- SuperLearner(Y = as.matrix(out_vec), X = out_pred_mat, newX = out_x0_mat, SL.library = aipw_lib)
#   # # pred_x0 <- c(learnOut_0$SL.predict)
#   # 
#   # psm <- SuperLearner(Y = as.matrix(out_vec), X = out_pred_mat, SL.library = aipw_lib)
#   # pred_x1 <- predict(psm , newdata = out_x1_mat)$pred
#   # pred_x0 <- predict(psm , newdata = out_x0_mat)$pred
#   # 
#   # gset <- gset %>% mutate(ps_u = 1-pred_ps, ps_t = pred_ps, pred_u = pred_x0, pred_t = pred_x1)
#   # 
#   # 
#   # # ps_u = (1-denom_fit$fitted.values), ps_t = denom_fit$fitted.values,
#   # # pred_u = predict(gmodel, data.frame(C1 = gset$C1, C2 = gset$C2, C3 = gset$C3,
#   # #                                     C4 = gset$C4, C5 = gset$C5, A = rep(0, nrow(gset)))),
#   # # pred_t = predict(gmodel, data.frame(C1 = gset$C1, C2 = gset$C2, C3 = gset$C3,
#   # #                                     C4 = gset$C4, C5 = gset$C5, A = rep(1, nrow(gset)))))
#   # 
#   # aipw <- RCAL::ate.aipw(y = gset$Y, tr = gset$A, mfp = cbind(gset$ps_u, gset$ps_t), mfo = cbind(gset$pred_u, gset$pred_t))
#   # ATE <- aipw$diff.est[2]
#   # SE <- sqrt(aipw$diff.var[2])
#   # TYPE <- "AIPW"
#   # res <- cbind(ATE, SE, TYPE)
#   
#   # 
#   # # DC-AIPW
#   # outcome <- as.character(as.formula(outForm)[[2]])
#   # exposure <- as.character(as.formula(expForm)[[2]])
#   # DCAIPW <- DCDR_Multiple(data=gset, exposure=exposure, outcome=outcome,
#   #               covarsT=all.vars(as.formula(expForm))[-1],
#   #               covarsO=all.vars(as.formula(outForm))[-1],
#   #               learners=aipw_lib,
#   #               control=SuperLearner.CV.control(V=5), num_cf=10)
#   # 
#   # ATE <- DCAIPW$rd
#   # SE <- DCAIPW$mvd
#   # TYPE <- "DC-AIPW"
#   # res <- cbind(ATE, SE, TYPE)
# 
#   
#   # IPW + GLM
#   num <- summary(glm(data = gset, formula = A ~ 1, family = "gaussian"))$coefficient[1,1]
#   denom_fit <- glm(data = gset,
#                    formula = as.formula(expForm),
#                    family = "binomial")
#   # gset <- select(gset, -denom)
#   gset <- gset %>% add_column(denom = denom_fit$fitted.values) %>%
#     mutate(wt = if_else(A == 1, num/denom, (1-num)/(1-denom))) #Stabilized IPTW
#   gset <- gset %>% mutate(wt2 = case_when(wt > quantile(gset$wt, 0.95) ~ quantile(gset$wt, 0.95),
#                                           wt < quantile(gset$wt, 0.05) ~ quantile(gset$wt, 0.05),
#                                           T ~ wt)) # Truncated IPTW
#   model_glm <- glm(data = gset, weight = wt2,
#                    formula = as.formula(outForm),
#                    family = "gaussian")
#   mean(gset[which(gset$A1==1),]$Y5)
#   mean(gset[which(gset$A1==0),]$Y5)
#   ATE <- summary(model_glm)$coefficients[2,1]
#   SE <- summary(model_glm)$coefficients[2,2]
#   TYPE <- "GLM + IPW"
#   res <- cbind(ATE, SE, TYPE)
# 
#   
#   
#   # CV-TMLE default specification
#   # nodes <- list(W = c("C1","C2","C3","C4","C5"), A = "A", Y = "Y")
#   # tmle3_autofit <- tmle3(tmle_ATE(treatment_level = 1, control_level = 0), data.table::copy(tset), nodes, learner_list)
#   # ATE <- tmle3_autofit$summary$tmle_est
#   # SE <- tmle3_autofit$summary$se
#   
#   
#   
#   # CV-TMLE full specification
#   tmle_task <- tmle3_Task$new(tset, npsem = npsem)
#   factor_list <- list(define_lf(LF_emp, "Z"), define_lf(LF_fit, "A", learner = tmle_lib), define_lf(LF_fit, "Y", learner = tmle_lib))
#   # likelihood_def <- Likelihood$new(factor_list)
#   # likelihood <- likelihood_def$train(tmle_task)
#   
#   learner_list <- list(A = tmle_lib, Y = tmle_lib)
#   likelihood <- ate_spec$make_initial_likelihood(
#     tmle_task,
#     learner_list
#   )
#   
#   ate_params <- list(Param_ATE$new(likelihood, define_lf(LF_static, "A", value = 1), define_lf(LF_static, "A", value = 0)))
#   
#   updater <- tmle3_Update$new()
#   targeted_likelihood <- Targeted_Likelihood$new(likelihood, updater)
#   
#   tmle3_fit <- fit_tmle3(tmle_task, targeted_likelihood, ate_params, updater)
#   ATE <- tmle3_fit$summary$tmle_est * (max(tset$Y)-min(tset$Y)) # backtransform ATE
#   SE <- tmle3_fit$summary$se * (max(tset$Y)-min(tset$Y)) # backtransform SE
# 
#   TYPE <- "CV-TMLE"
# 
#   res <- cbind(ATE, SE, TYPE)
#   # res <- rbind(res, cbind(ATE, SE, TYPE))
# 
#   # Manual TMLE
#   outcome <- as.character(as.formula(outForm)[[2]])
#   exposure <- as.character(as.formula(expForm)[[2]])
#   TMLE <- TMLE(data=tset, exposure=exposure, outcome=outcome,
#                           covarsT=all.vars(as.formula(expForm))[-1],
#                           covarsO=all.vars(as.formula(outForm))[-1],
#                           learners=aipw_lib,
#                           control=SuperLearner.CV.control(V=5))
# 
#   ATE <- TMLE$rd
#   SE <- sqrt(TMLE$vd)
#   TYPE <- "TMLE"
#   # res <- cbind(ATE, SE, TYPE)
#   res <- rbind(res, cbind(ATE, SE, TYPE))
#   
#   
#   # # DC-TMLE
#   # outcome <- as.character(as.formula(outForm)[[2]])
#   # exposure <- as.character(as.formula(expForm)[[2]])
#   # DCTMLE <- DCTMLE_Multiple(data=tset, exposure=exposure, outcome=outcome,
#   #                         covarsT=all.vars(as.formula(expForm))[-1],
#   #                         covarsO=all.vars(as.formula(outForm))[-1],
#   #                         learners=aipw_lib,
#   #                         control=SuperLearner.CV.control(V=5), num_cf=10)
#   # 
#   # ATE <- DCTMLE$rd
#   # SE <- sqrt(DCTMLE$mvd)
#   # TYPE <- "DC-TMLE"
#   # res <- cbind(ATE, SE, TYPE)
#   return(res)
# }


# meta_res <- NULL # intialize the results table
source("20200705-DCDR-Functions.R")
source("20200720-Algos-code.R")

##################################
## PARALLELIZE ANALYSES
##################################
{
set.seed(42782)
tic()
#N_sims <- 1000
N_sims <- 10
#ss <- nrow(tmp3_72)
ss <- 600
# specify which set of learners for SL
tmle_lib <- lrnr_SL
#tmle_lib <- lrnr_SL_param
aipw_lib <- SL.lib
#aipw_lib <- SL.param

# regression models for GLM / AIPW
expForm <- "A ~ C1 + C2 + C3 + C4 + C5"
outForm <- "Y ~ A + C1 + C2 + C3 + C4 + C5"

boot1 <- foreach(i = 1:N_sims) %dopar% {
  require(tidyverse)
  require(tmle3)
  require(sl3)
  require(SuperLearner)
  
  # Initialize dataset
  set1 <- as_tibble(cbind(C1 = sim_boots2[[i]]$C1[1:ss],
                    C2 = sim_boots2[[i]]$C2[1:ss],
                    C3 = sim_boots2[[i]]$C3[1:ss],
                    C4 = sim_boots2[[i]]$C4[1:ss],
                    C5 = sim_boots2[[i]]$C5[1:ss],
                    A = sim_boots2[[i]]$A[1:ss],
                    Y = sim_boots2[[i]]$Y[1:ss]))
  tset <- set1 %>% mutate(YT = (Y-min(set1$Y))/(max(set1$Y)- min(set1$Y))) # generate a bounded Y for TMLE
  getRES(set1, tset, aipw_lib, tmle_lib, short_tmle_lib,
         doAIPW=0, doDCAIPW=0,
         doIPW = 0, 
         doTMLE=1, doManuTMLE=1, doShortTMLE = 0,
         doDCTMLE=0
  )
}
toc()
}

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

sim_res
out_path <- "/Users/garethalex/Desktop/HuangGroup/cvtmle_plasmode/Data/"
boot1.out <- data.frame(matrix(unlist(tmp),ncol=3))
fac.to.num <- function(f) as.numeric(levels(f))[f]
boot1.out[,1:2] <- lapply(boot1.out[,1:2], fac.to.num)
boot1.out %>% group_by(X3)
sim_res
# writexl::write_xlsx(boot1.out, paste0(out_path,"sim_results_for_DCDR_detail.xlsx"))
}







{
pos <- max(sim_corr1$ub)
sim_corr1 %>% 
  ggplot(aes(x = iter/2, y = ATE, color = TYPE)) +
  geom_point() + geom_errorbar(aes(ymin = lb, ymax = ub)) + 
  geom_hline(aes(yintercept = Effect_Size)) + 
  geom_hline(data = filter(sim_corr1, substring(TYPE,1,3) %in% c("AIP", "GLM")), aes(yintercept = mean(ATE)), linetype = "dashed") +
  geom_hline(data = filter(sim_corr1, substring(TYPE,1,2) == "CV"), aes(yintercept = mean(ATE)), linetype = "dotted") +
  geom_text(data = filter(sim_corr1, substring(TYPE,1,3) %in% c("AIP", "GLM")),
            aes(N_sims, pos, hjust = "right",
                label = paste("\n\nMean bias = ", as.character(round(mean(bias), 3)),
                              "\nMedian bias = ", as.character(round(median(bias), 3)),
                              "\nMSE = ", as.character(round(var(ATE) + mean(bias)^2, 3)),
                              "\nCoverage = ", as.character(round(sum(lb <= Effect_Size & ub >= Effect_Size)/N_sims, 3)) ))) +
  geom_text(data = filter(sim_corr1, substring(TYPE,1,2) == "CV"),
            aes(N_sims, pos, hjust = "right",
                label = paste("\n\nMean bias = ", as.character(round(mean(bias), 3)),
                              "\nMedian bias = ", as.character(round(median(bias), 3)),
                              "\nMSE = ", as.character(round(var(ATE) + mean(bias)^2, 3)),
                              "\nCoverage = ", as.character(round(sum(lb <= Effect_Size & ub >= Effect_Size)/N_sims, 3)) ))) +
  labs(x = "iteration", color = "Estimator") +
  facet_wrap(~TYPE) +
  theme(legend.position = "none")

# append results
meta_res <- rbind(meta_res, 
                  cbind(N = ss, iter = N_sims, 
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

out_path <- "/Users/garethalex/Desktop/HuangGroup/cvtmle_plasmode/Data/"
writexl::write_xlsx(data.frame(meta_res), paste0(out_path,"sim_results_for_DCDR.xlsx"))
boot1.out <- data.frame(matrix(unlist(boot1), nrow=length(boot1), byrow=T))
fac.to.num <- function(f) as.numeric(levels(f))[f]
boot1.out[,1:2] <- lapply(boot1.out[,1:2], fac.to.num)
boot1.out
# writexl::write_xlsx(boot1.out, paste0(out_path,"sim_results_for_DCDR_detail.xlsx"))


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
