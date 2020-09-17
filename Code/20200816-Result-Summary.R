
# run.est <- function(){
#   set.seed(42782)
#   tic()
#   N_sims <- 10# this should <= plas_sim_N
#   # DC implementation
#   source("20200705-DCDR-Functions.R")
#   # getRES function is now relocated to a separate file (20200720-Algos-code.R)
#   source("20200720-Algos-code.R")
#   
#   # regression models for GLM / AIPW
#   if (estimateWithMore == T){
#     # out_path <- "/Users/garethalex/Desktop/HuangGroup/cvtmle_plasmode/Data/"
#     plas_org <- haven::read_dta(paste0("./plas_data.dta"))
#     vars <- names(plas_org[3:333])[1:p.est]
#     p <- p.est
#   }
#   reg.formulas <- make.formula("Y", "A",ver = data.ver,sims.ver = sims.ver,vars=vars,p=p.est)
#   expForm <- reg.formulas$expForm
#   outForm <- reg.formulas$outForm
#   
#   # specify which set of learners for SL
#   # ### NON-SMOOTH
#   # short_tmle_lib <- SL_param_list
#   # tmle_lib <- lrnr_SL
#   # aipw_lib <- SL.lib
#   # 
#   # # SMOOTH
#   # short_tmle_lib <- SL_list
#   # tmle_lib <- lrnr_SL_param
#   # aipw_lib <- SL.param
#   
#   
#   # set .errorhandling="remove" if want to discard
#   # set .errorhandling="stop" by default
#   
#   # boot1 <- foreach(i = 1:N_sims,.errorhandling="stop") %dopar% {
#   boot1 <- foreach(i = 1:N_sims,.errorhandling="remove") %dopar% {
#     require(tidyverse)
#     require(tmle3)
#     require(sl3)
#     require(SuperLearner)
#     
#     # Initialize dataset
#     if (sims.ver == "plas" | sims.ver =="5var.then.plas"){
#       # i <- 11
#       plas_data <- cbind(id = plas_sims$Sim_Data[i],
#                          A = plas_sims$Sim_Data[i + (2*plas_sim_N)],
#                          Y = plas_sims$Sim_Data[i + plas_sim_N])
#       colnames(plas_data) <- c("id", "A", "Y")
#       set1 <- suppressMessages(left_join(as_tibble(plas_data), as_tibble(plas))) #dplyr::select(as_tibble(plas), -Y5, -A1))) # add covariates
#       tset <- set1 %>% mutate(YT = (Y-min(set1$Y))/(max(set1$Y)- min(set1$Y)))
#       
#     }else{
#       # Initialize dataset
#       # i<-2
#       ss <- 600
#       set1 <- as_tibble(cbind(C1 = sim_boots[[i]]$C1[1:ss],
#                               C2 = sim_boots[[i]]$C2[1:ss],
#                               C3 = sim_boots[[i]]$C3[1:ss],
#                               C4 = sim_boots[[i]]$C4[1:ss],
#                               C5 = sim_boots[[i]]$C5[1:ss],
#                               A = sim_boots[[i]]$A[1:ss],
#                               Y = sim_boots[[i]]$Y[1:ss]))
#       tset <- set1 %>% mutate(YT = (Y-min(set1$Y))/(max(set1$Y)- min(set1$Y))) # generate a bounded Y for TMLE
#     }
#     getRES(set1, tset, aipw_lib, tmle_lib, short_tmle_lib,
#            doIPW = 0, doLASSO=0,
#            doAIPW=1, doDCAIPW=0,
#            doManuTMLE=0, doShortTMLE = 0,
#            doDCTMLE=0,
#            num_cf=5,
#            #control=list()
#            control=SuperLearner.CV.control(V=2)
#     )
#   }
#   toc()
# }

Effect_Size <- 6.6
summarise.res <- function(boot1){
  tmp <- NULL
  N_boot <- length(boot1)
  for(i in 1:N_boot){
    tmp <- rbind(tmp, boot1[[i]])
  }
  
  sim_corr1 <- as_tibble(tmp) %>% 
    mutate(ATE = as.double(ATE), SE = as.double(SE),
           lb = as.double(ATE) - 1.96*as.double(SE), 
           ub = as.double(ATE) + 1.96*as.double(SE)) %>% 
    add_column(iter = c(1:nrow(tmp))) %>% mutate(bias = as.double(ATE) - Effect_Size)
  
  #write_csv(plas_corr1, paste0(out_path, Sys.Date(), "-plasmode_sim_CV-TMLE_GLM_IPW_1000.csv"))
  
  sim_res <- sim_corr1 %>% group_by(TYPE) %>% 
    summarize(b=length(ATE), mu_ATE = mean(ATE), med_ATE = median(ATE), 
              mu_SE = mean(SE), med_SE = median(SE), 
              mu_bias = mean(bias), med_bias = median(bias), 
              var = var(ATE), MSE = var + mu_bias^2,
              coverage = sum(lb <= Effect_Size & ub >= Effect_Size)/N_boot)
  
  # out_path <- "/Users/garethalex/Desktop/HuangGroup/cvtmle_plasmode/Data/"
  # boot1.out <- data.frame(matrix(unlist(tmp),ncol=3))
  # fac.to.num <- function(f) as.numeric(levels(f))[f]
  # boot1.out[,1:2] <- lapply(boot1.out[,1:2], fac.to.num)
  # boot1.out %>% group_by(X3)
  # pos <- max(sim_corr1$ub) 
  
  
  # Visualize
  sim_plot <- sim_corr1 %>% 
    ggplot(aes(x = iter, y = ATE, color = TYPE)) +
    geom_point() + geom_errorbar(aes(ymin = lb, ymax = ub)) + 
    geom_hline(aes(yintercept = Effect_Size)) + 
    # geom_hline(data = filter(sim_corr1, substring(TYPE,1,3) %in% c("AIP", "GLM")), aes(yintercept = mean(ATE)), linetype = "dashed") +
    # geom_hline(data = filter(sim_corr1, substring(TYPE,1,2) == "CV"), aes(yintercept = mean(ATE)), linetype = "dotted") +
    geom_hline(data = sim_corr1, aes(yintercept = mean(ATE)), linetype = "dotted") +
    # geom_text(data = filter(sim_corr1, substring(TYPE,1,3) %in% c("AIP", "GLM")),
    #           aes(N_boot, pos, hjust = "right",
    #               label = paste("\n\nMean bias = ", as.character(round(mean(bias), 3)),
    #                             "\nMedian bias = ", as.character(round(median(bias), 3)),
    #                             "\nMSE = ", as.character(round(var(ATE) + mean(bias)^2, 3)),
    #                             "\nCoverage = ", as.character(round(sum(lb <= Effect_Size & ub >= Effect_Size)/N_boot, 3)) ))) +
    # geom_text(data = filter(sim_corr1, substring(TYPE,1,2) == "CV"),
    #           aes(N_boot, pos, hjust = "right",
    #               label = paste("\n\nMean bias = ", as.character(round(mean(bias), 3)),
    #                             "\nMedian bias = ", as.character(round(median(bias), 3)),
    #                             "\nMSE = ", as.character(round(var(ATE) + mean(bias)^2, 3)),
  #                             "\nCoverage = ", as.character(round(sum(lb <= Effect_Size & ub >= Effect_Size)/N_boot, 3)) ))) +
  labs(x = "iteration", color = "Estimator") +
    facet_wrap(~TYPE) +
    theme(legend.position = "none")
  # sim_plot
  return(list(sim_res,sim_plot))
}
