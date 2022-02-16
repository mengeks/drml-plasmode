library(ggplot2)
library(tidyverse)
# setwd("~/Desktop/HuangGroup/cvtmle_plasmode/Code")
# setwd("~/Desktop/HuangGroup/cvtmle_plasmode/Code")
Effect_Size <- 6.6
# Effect_Size <- 5.51788
# Effect_Size <- 1.190723
source("./20200816-Result-Summary.R")
is.timing <- T; tim.str <- ifelse(is.timing,"-timing","")


for (i in 1:length(boot1)){
  colnames(boot1[[i]]) <- c("ATE", "SE", "TYPE", "t", "no_cores")
}

# give.summary.res(res.path, "med_bias")
summarise.res(boot1, Effect_Size=Effect_Size)


tmp <- NULL
N_boot <- length(boot1)
for(i in 1:N_boot){
  tmp <- rbind(tmp, boot1[[i]])
}
if ("t" %in% colnames(tmp)){
  t.vec <- as.double(as_tibble(tmp)$t)
}else{
  t.vec <- NA
}
sim_corr1 <- as_tibble(tmp) %>% 
  mutate(ATE = as.double(ATE), SE = as.double(SE),
         lb = as.double(ATE) - 1.96*as.double(SE), 
         ub = as.double(ATE) + 1.96*as.double(SE),
         t = t.vec, kl=as.double(no_cores)) %>% 
  add_column(iter = c(1:nrow(tmp))) %>% mutate(bias = as.double(ATE) - Effect_Size)

#write_csv(plas_corr1, paste0(out_path, Sys.Date(), "-plasmode_sim_CV-TMLE_GLM_IPW_1000.csv"))
plot(sim_corr1$kl,abs(sim_corr1$bias))
plot(sim_corr1$kl,abs(sim_corr1$bias), ylim=c(0,1))

sim_res <- sim_corr1 %>% group_by(TYPE) %>% 
  summarize(b=length(ATE), med_t=median(t), mu_ATE = mean(ATE), med_ATE = median(ATE), 
            mu_SE = mean(SE), med_SE = median(SE), 
            mu_bias = mean(bias), med_bias = median(bias), 
            var = var(ATE), MSE = var + mu_bias^2,
            coverage = sum(lb <= Effect_Size & ub >= Effect_Size)/N_boot)