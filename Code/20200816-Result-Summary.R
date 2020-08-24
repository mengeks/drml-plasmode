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
    summarize(mu_ATE = mean(ATE), med_ATE = median(ATE), 
              mu_SE = mean(SE), med_SE = median(SE), 
              mu_bias = mean(bias), med_bias = median(bias), 
              var = var(ATE), MSE = var + mu_bias^2,
              coverage = sum(lb <= Effect_Size & ub >= Effect_Size)/N_boot)
  
  out_path <- "/Users/garethalex/Desktop/HuangGroup/cvtmle_plasmode/Data/"
  boot1.out <- data.frame(matrix(unlist(tmp),ncol=3))
  fac.to.num <- function(f) as.numeric(levels(f))[f]
  boot1.out[,1:2] <- lapply(boot1.out[,1:2], fac.to.num)
  boot1.out %>% group_by(X3)
  pos <- max(sim_corr1$ub) 
  
  
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
