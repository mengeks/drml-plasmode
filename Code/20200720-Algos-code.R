

#expForm <- "A ~ L0.a + L0.b + L0.c + L0.d + L0.e + L0.f + L0.g + L0.h + L0.i + L0.j + L0.k"
#outForm <- "Y ~ A + L0.a + L0.b + L0.c + L0.d + L0.e + L0.f + L0.g + L0.h + L0.i + L0.j + L0.k"

###############################
# FUNCTION TO RUN AIPW, CV-TMLE, and IPW+GLM
###############################
getRES <- function(gset, tset, aipw_lib = NULL, tmle_lib = NULL, short_tmle_lib=NULL,
                   doAIPW =0, doDCAIPW=0,
                   doIPW = 0,
                   doTMLE = 0, doManuTMLE=0, doShortTMLE = 0,
                   doDCTMLE=0,
                   num_cf=5,
                   control){
  res <- NULL
  
  # AIPW
  if (doAIPW == 1){
    # AIPW
    outcome <- as.character(as.formula(outForm)[[2]])
    exposure <- as.character(as.formula(expForm)[[2]])
    out_vec <- pull(select(gset, outcome))
    exp_vec <- pull(select(gset, exposure))
    
    # denom_fit <- glm(data = gset, formula = as.formula(expForm), family = "binomial",fit=FALSE)
    # gmodel <- glm(data = gset, formula = as.formula(outForm), family = "gaussian",fit=FALSE)
    #
    # expVars <- names(denom_fit$coefficients)[2:length(denom_fit$coefficients)]
    # exp_pred_mat <- data.frame(select(gset, expVars))
    exp_pred_mat <- model.matrix(object = as.formula(expForm), data=gset)[,-1]
    exp_pred_mat <- data.frame(exp_pred_mat)
    suppressMessages(learnPS <- SuperLearner(Y = exp_vec, X = exp_pred_mat, SL.library = aipw_lib, family = binomial()))
    pred_ps <- c(learnPS$SL.predict)
    
    # outVars <- names(gmodel$coefficients)[2:length(gmodel$coefficients)]
    # out_pred_mat <- data.frame(select(gset, outVars))
    out_pred_mat <- model.matrix(object = as.formula(outForm), data=gset)[,-1]
    out_pred_mat <- data.frame(out_pred_mat)
    out_x1_mat <- out_x0_mat <- out_pred_mat
    out_x1_mat[exposure] <- rep(1, length(exp_vec))
    out_x0_mat[exposure] <- rep(0, length(exp_vec))
    
    # learnOut_1 <- SuperLearner(Y = as.matrix(out_vec), X = out_pred_mat, newX = out_x1_mat, SL.library = aipw_lib)
    # pred_x1 <- c(learnOut_1$SL.predict)
    # learnOut_0 <- SuperLearner(Y = as.matrix(out_vec), X = out_pred_mat, newX = out_x0_mat, SL.library = aipw_lib)
    # pred_x0 <- c(learnOut_0$SL.predict)
    
    psm <- SuperLearner(Y = as.matrix(out_vec), X = out_pred_mat, SL.library = aipw_lib)
    pred_x1 <- predict(psm , newdata = out_x1_mat)$pred
    pred_x0 <- predict(psm , newdata = out_x0_mat)$pred
    
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
  
  if (doDCAIPW==1){
    
    
    # DC-AIPW
    outcome <- as.character(as.formula(outForm)[[2]])
    exposure <- as.character(as.formula(expForm)[[2]])
    DCAIPW <- DCDR_Multiple(data=gset, exposure=exposure, outcome=outcome,
                            covarsT=all.vars(as.formula(expForm))[-1],
                            covarsO=all.vars(as.formula(outForm))[-1],
                            learners=aipw_lib,
                            control=control, 
                            num_cf=num_cf)
    
    ATE <- DCAIPW$rd
    SE <- DCAIPW$mvd
    TYPE <- "DC-AIPW"
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
  
  if (doShortTMLE == 1){
    node_list <- list(W = Zvars,
                      A = "A",
                      Y = "YT")
    
    ate_spec <- tmle_ATE(
      treatment_level = 1,
      control_level = 0
    )
    
    learner_list <- short_tmle_lib
    tmle_fit <- tmle3(ate_spec,tset, node_list, learner_list)
    ATE <- tmle_fit$summary$psi_transformed * (max(tset$Y)-min(tset$Y))
    SE <- tmle_fit$summary$se * (max(tset$Y)-min(tset$Y))
    TYPE <- "Short-TMLE"
    res <- rbind(res, cbind(ATE, SE, TYPE))
  }
  
  if (doManuTMLE==1){
    # Manual TMLE
    outcome <- as.character(as.formula(outForm)[[2]])
    exposure <- as.character(as.formula(expForm)[[2]])
    TMLE <- TMLE(data=tset, exposure=exposure, outcome=outcome,
                 covarsT=all.vars(as.formula(expForm))[-1],
                 covarsO=all.vars(as.formula(outForm))[-1],
                 learners=aipw_lib,
                 control=control)
    
    ATE <- TMLE$rd
    SE <- sqrt(TMLE$vd)
    TYPE <- "ManuTMLE"
    res <- rbind(res, cbind(ATE, SE, TYPE))
  }
  
  if (doDCTMLE==1){
    # DC-TMLE
    outcome <- as.character(as.formula(outForm)[[2]])
    exposure <- as.character(as.formula(expForm)[[2]])
    DCTMLE <- DCTMLE_Multiple(data=tset, exposure=exposure, outcome=outcome,
                              covarsT=all.vars(as.formula(expForm))[-1],
                              covarsO=all.vars(as.formula(outForm))[-1],
                              learners=aipw_lib,
                              control=control, num_cf=num_cf)
    
    ATE <- DCTMLE$rd
    SE <- sqrt(DCTMLE$mvd)
    TYPE <- "DC-TMLE"
    res <- rbind(res, cbind(ATE, SE, TYPE))
  }
  
  return(res)
}