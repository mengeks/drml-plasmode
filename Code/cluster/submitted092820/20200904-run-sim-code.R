{
  set.seed(42782)
  tic()
  
  # DC implementation
  source("20200705-DCDR-Functions.R")
  # getRES function is now relocated to a separate file (20200720-Algos-code.R)
  source("20200720-Algos-code.R")
  
  # regression models for GLM / AIPW
  if (estimateWithMore == T){
    out_path <- "/Users/garethalex/Desktop/HuangGroup/cvtmle_plasmode/Data/"
    plas_org <- haven::read_dta(paste0(out_path,"plas_data.dta"))
    vars <- names(plas_org[3:333])[1:p.est]
    p <- p.est
  }
  reg.formulas <- make.formula("Y", "A",ver = data.ver,sims.ver = sims.ver,vars=vars,p=p.est)
  expForm <- reg.formulas$expForm
  outForm <- reg.formulas$outForm
  
  
  
  
  # set .errorhandling="remove" if want to discard
  # set .errorhandling="stop" by default
  
  # boot1 <- foreach(i = 1:N_sims,.errorhandling="stop") %dopar% {
  boot1 <- foreach(i = 1:N_sims,.errorhandling=errorhandling) %dopar% {
    require(tidyverse)
    require(tmle3)
    require(sl3)
    require(SuperLearner)
    
    # Initialize dataset
    if (sims.ver == "plas" | sims.ver =="5var.then.plas"){
      # i <- 11
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
           doIPW = doIPW, doLASSO=doLASSO,
           doAIPW=doAIPW, 
           doDCAIPW=doDCAIPW,
           doManuTMLE=doManuTMLE, doShortTMLE = doShortTMLE,
           doDCTMLE=doDCTMLE,
           num_cf=num_cf,
           #control=list()
           control=control
    )
  }
  toc()
}