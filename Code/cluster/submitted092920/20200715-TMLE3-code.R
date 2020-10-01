library(data.table)
library(tidyverse)
library(tmle3)
library(sl3)
washb_data <- fread("https://raw.githubusercontent.com/tlverse/tlverse-data/master/wash-benefits/washb_data.csv",
                    stringsAsFactors = TRUE)
head(washb_data)
washb_data$tr

# Initialize dataset
i<-1
plas_data <- cbind(id = plas_sims$Sim_Data[i],
                   A = plas_sims$Sim_Data[i + (2*plas_sim_N)],
                   Y = plas_sims$Sim_Data[i + plas_sim_N])
colnames(plas_data) <- c("id", "A", "Y")
set1 <- suppressMessages(left_join(as_tibble(plas_data), as_tibble(plas))) #dplyr::select(as_tibble(plas), -Y5, -A1))) # add covariates
# mean(set1[which(set1$A1==1),]$Y) - mean(set1[which(set1$A1==0),]$Y)
tset <- set1 %>% mutate(YT = (Y-min(set1$Y))/(max(set1$Y)- min(set1$Y)))

# node_list <- list(W = c("C1", "C2", "C3", "C4", "C5"),
#               A = "A",
#               Y = "YT")
{
  tic()
node_list <- list(W = Zvars,
              A = "A",
              Y = "YT")

ate_spec <- tmle_ATE(
  treatment_level = 1,
  control_level = 0
)

learner_list <- list(A = lrnr_SL, Y = lrnr_SL)
tmle_fit <- tmle3(ate_spec,tset, node_list, learner_list)
ATE <- tmle_fit$summary$psi_transformed * (max(tset$Y)-min(tset$Y))
SE <- tmle_fit$summary$se * (max(tset$Y)-min(tset$Y))
toc()
}





tmle_task <- ate_spec$make_tmle_task(tset, node_list)
tmle_task$npsem
initial_likelihood <- ate_spec$make_initial_likelihood(
  tmle_task,
  learner_list
)
print(initial_likelihood)
initial_likelihood$get_likelihoods(tmle_task)
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood)
targeted_likelihood_no_cv <-
  Targeted_Likelihood$new(initial_likelihood,
                          updater = list(cvtmle = FALSE)
  )
tmle_params <- ate_spec$make_params(tmle_task, targeted_likelihood)
# print(tmle_params)
tmle_fit_manual <- fit_tmle3(
  tmle_task, targeted_likelihood, tmle_params,
  targeted_likelihood$updater
)
print(tmle_fit_manual)
estimates <- tmle_fit_manual$summary$tmle_est * (max(tset$Y)-min(tset$Y))
print(estimates)

