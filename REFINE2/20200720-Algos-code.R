
expit <- function(x){
  1/(1+exp(-x))
}

#expForm <- "A ~ L0.a + L0.b + L0.c + L0.d + L0.e + L0.f + L0.g + L0.h + L0.i + L0.j + L0.k"
#outForm <- "Y ~ A + L0.a + L0.b + L0.c + L0.d + L0.e + L0.f + L0.g + L0.h + L0.i + L0.j + L0.k"
##########################################
# Initialize necessary parameters for estimator
##########################################

{
  # TMLE parameters
  #SL.lib <- c("SL.randomForest", "SL.xgboost", "SL.nnet", "SL.glm", "SL.glmnet", "SL.polymars")
  #SL.lib <- list("SL.randomForest", "SL.xgboost", "SL.nnet", "SL.glm", c("SL.glmnet", "All"), c("SL.polymars", "All"))
  SL.lib <- list("SL.randomForest", "SL.xgboost", "SL.glm", c("SL.polymars", "All"))
  SL.lib.tmle <- c("SL.randomForest", "SL.xgboost", "SL.glm", "SL.polymars")
  SL.param <- c("SL.glm", "SL.glmnet", "SL.polymars")
  SL.lib.both <- c(SL.lib.tmle,SL.param)
  
  # Specify the NP-SEM for the TMLE - including bounded, tranformed Y ("YT")
  # Zvars <- ifelse(data.ver=="FULL",
  #                 vars,
  #                 c("L0.a", "L0.b", "L0.c", "L0.d", "L0.e", "L0.f", "L0.g", "L0.h", "L0.i", "L0.j", "L0.k"))
  # npsem <- list(define_node("Z", Zvars),
  #               #c("L0.a", "L0.b", "L0.c", "L0.d", "L0.e", "L0.f", "L0.g", "L0.h", "L0.i", "L0.j", "L0.k")),
  #               define_node("A", c("A"), c("Z")),
  #               define_node("Y", c("YT"), c("A", "Z")))
  
  # Specify the learners for CV-TMLE
  # lrnr_SL <- make_learner(Lrnr_pkg_SuperLearner, SL.lib.tmle)
  # lrnr_SL_param <- make_learner(Lrnr_pkg_SuperLearner, SL.param)
  # lrnr_SL_both <- make_learner(Lrnr_pkg_SuperLearner, SL.lib.both)
  
  # lrnr_glm <- make_learner(Lrnr_glm_fast)
  # lrnr_mean <- make_learner(Lrnr_mean)
  
  # SL_list <- list(Y = lrnr_SL, A = lrnr_SL)
  # SL_param_list <- list(Y = lrnr_SL_param, A = lrnr_SL_param)
  # glm_list <- list(Y = lrnr_glm, A = lrnr_glm)
  
  # short_tmle_lib <- SL_param_list
  # tmle_lib <- lrnr_SL
  aipw_lib <- SL.lib
}

###############################
# FUNCTION TO RUN AIPW, CV-TMLE, and IPW+GLM
###############################
getRES <- function(gset, tset, aipw_lib = NULL, tmle_lib = NULL, short_tmle_lib=NULL,
                   doAIPW =0, doDCAIPW=0,
                   doIPW = 0, doLASSO=0,
                   doTMLE = 0, doManuTMLE=0, doShortTMLE = 0,
                   doDCTMLE=0,doGComp=0,
                   doManuTMLE_Auto=0, doAIPW_Auto=0,
                   doDRTMLE=0,
                   doAIPW_Ensemble=0,
                   num_cf=5,
                   control,parallel=T,
                   expForm, outForm,no_cores=8){
  res <- NULL
  
  # AIPW
  if (doAIPW == 1){
    # AIPW
    ptm0 <- proc.time()
    outcome <- as.character(as.formula(outForm)[[2]])
    exposure <- as.character(as.formula(expForm)[[2]])
    out_vec <- pull(select(gset, outcome))
    exp_vec <- pull(select(gset, exposure))
    
    exp_pred_mat <- model.matrix(object = as.formula(expForm), data=gset)[,-1]
    exp_pred_mat <- data.frame(exp_pred_mat)
    # suppressMessages(learnPS <- SuperLearner(Y = exp_vec, X = exp_pred_mat, 
    #                                          SL.library = c(aipw_lib,rf_names), family = binomial(),
    #                                          cvControl=control))
    learnPS <- SuperLearner(Y = exp_vec, X = exp_pred_mat, 
                            SL.library = aipw_lib, family = binomial(),
                            cvControl=control)
    
    # observe({print(head(exp_pred_mat))})
    # SuperLearner(Y=as.matrix(data[, exposure])[,1], X=data.frame(data[, covarsT]), 
    #              family=binomial(), SL.library=learners, cvControl=control)
    pred_ps <- c(learnPS$SL.predict)
    # ptmps <- proc.time()
    # print(paste0("pscore fitting takes ", round((ptmps - ptm0)[3],3), "s"))

    out_pred_mat <- model.matrix(object = as.formula(outForm), data=gset)[,-1]
    out_pred_mat <- data.frame(out_pred_mat)
    out_x1_mat <- out_x0_mat <- out_pred_mat
    out_x1_mat[exposure] <- rep(1, length(exp_vec))
    out_x0_mat[exposure] <- rep(0, length(exp_vec))
    
    psm <- SuperLearner(Y = as.matrix(out_vec), X = out_pred_mat, SL.library = aipw_lib,cvControl=control)
    # ptmor <- proc.time()
    # print(paste0("OR fitting takes ", round((ptmor - ptmps)[3],3), "s"))
    pred_x1 <- predict(psm , newdata = out_x1_mat)$pred
    pred_x0 <- predict(psm , newdata = out_x0_mat)$pred
    
    gset <- gset %>% mutate(ps_u = 1-pred_ps, ps_t = pred_ps, pred_u = pred_x0, pred_t = pred_x1)
    
    
    aipw <- RCAL::ate.aipw(y = gset$Y, tr = gset$A, mfp = cbind(gset$ps_u, gset$ps_t), mfo = cbind(gset$pred_u, gset$pred_t))
    # ptmrcal <- proc.time()
    # print(paste0("RCAL takes ", round((ptmrcal - ptmor)[3],3), "s"))
    ATE <- aipw$diff.est[2]
    SE <- sqrt(aipw$diff.var[2])
    TYPE <- "AIPW"
    ptm1 <- proc.time()
    t <- round((ptm1 - ptm0)[3],3)
    res <- rbind(res, cbind(ATE, SE, TYPE,t, no_cores))
  }
  
  if (doAIPW_Auto == 1){
    # AIPW
    ptm0 <- proc.time()
    list_nTrees = c(50, 100, 1000) # try 5 different values for nTree
    list_mTry = c(10,20) # try 6 different values for mTry (default is sqrt(M), where M is the total number of X vars)
    rf_list <- list(); rf_names <- c(); idx <- 1
    for (i in 1:length(list_nTrees)){
      for (j in 1:length(list_mTry)){
        ntree <- list_nTrees[i]; mtry <- list_mTry[j]
        rf_list <- c(rf_list,function(...){SL.randomForest(..., ntree=ntree, mtry=mtry)})
        rf_names <- c(rf_names, paste0("SL.randomForest_ntree",ntree,"mtry",mtry))
        idx <- idx+1
      }
    }
    names(rf_list) <- rf_names
    assign("rf_list",rf_list, envir = .GlobalEnv)
    list2env(rf_list,globalenv())
    
    outcome <- as.character(as.formula(outForm)[[2]])
    exposure <- as.character(as.formula(expForm)[[2]])
    out_vec <- pull(select(gset, outcome))
    exp_vec <- pull(select(gset, exposure))
    
    exp_pred_mat <- model.matrix(object = as.formula(expForm), data=gset)[,-1]
    exp_pred_mat <- data.frame(exp_pred_mat)
    # suppressMessages(learnPS <- SuperLearner(Y = exp_vec, X = exp_pred_mat, 
    #                                          SL.library = c(aipw_lib,rf_names), family = binomial(),
    #                                          cvControl=control))
    learnPS <- SuperLearner(Y = exp_vec, X = exp_pred_mat, 
                                             SL.library = c(aipw_lib,rf_names), family = binomial(),
                                             cvControl=control)
    
    observe({print(head(exp_pred_mat))})
    pred_ps <- c(learnPS$SL.predict)
    # ptmps <- proc.time()
    # print(paste0("pscore fitting takes ", round((ptmps - ptm0)[3],3), "s"))
    
    out_pred_mat <- model.matrix(object = as.formula(outForm), data=gset)[,-1]
    out_pred_mat <- data.frame(out_pred_mat)
    out_x1_mat <- out_x0_mat <- out_pred_mat
    out_x1_mat[exposure] <- rep(1, length(exp_vec))
    out_x0_mat[exposure] <- rep(0, length(exp_vec))
    
    psm <- SuperLearner(Y = as.matrix(out_vec), X = out_pred_mat, 
                        SL.library = c(aipw_lib,rf_names),cvControl=control)
    # ptmor <- proc.time()
    # print(paste0("OR fitting takes ", round((ptmor - ptmps)[3],3), "s"))
    pred_x1 <- predict(psm , newdata = out_x1_mat)$pred
    pred_x0 <- predict(psm , newdata = out_x0_mat)$pred
    
    gset <- gset %>% mutate(ps_u = 1-pred_ps, ps_t = pred_ps, pred_u = pred_x0, pred_t = pred_x1)
    
    
    aipw <- RCAL::ate.aipw(y = gset$Y, tr = gset$A, mfp = cbind(gset$ps_u, gset$ps_t), mfo = cbind(gset$pred_u, gset$pred_t))
    # ptmrcal <- proc.time()
    # print(paste0("RCAL takes ", round((ptmrcal - ptmor)[3],3), "s"))
    ATE <- aipw$diff.est[2]
    SE <- sqrt(aipw$diff.var[2])
    TYPE <- "AIPW"
    ptm1 <- proc.time()
    t <- round((ptm1 - ptm0)[3],3)
    res <- rbind(res, cbind(ATE, SE, TYPE,t, no_cores))
  }
  
  if (doDCAIPW==1){
    
    ptm0 <- proc.time()
    # DC-AIPW
    outcome <- as.character(as.formula(outForm)[[2]])
    exposure <- as.character(as.formula(expForm)[[2]])
    DCAIPW <- DCDR_Multiple(data=gset, exposure=exposure, outcome=outcome,
                            covarsT=all.vars(as.formula(expForm))[-1],
                            covarsO=all.vars(as.formula(outForm))[-1],
                            learners=aipw_lib,
                            control=control, 
                            num_cf=num_cf,parallel=parallel)
    
    ATE <- DCAIPW$rd
    SE <- sqrt(DCAIPW$mvd)
    TYPE <- "DC-AIPW"
    ptm1 <- proc.time()
    t <- round((ptm1 - ptm0)[3],3)
    res <- rbind(res, cbind(ATE, SE, TYPE,t, no_cores))
  }
  
  # IPW + GLM
  if (doIPW == 1){
    # gset <- set1
    ptm0 <- proc.time()
    num <- summary(glm(data = gset, formula = A ~ 1, family = "gaussian"))$coefficient[1,1]
    denom_fit <- glm(data = gset,
                     formula = as.formula(expForm),
                     family = "binomial",
                     maxit=100)
    summary(denom_fit)
    gset <- gset %>% add_column(denom = denom_fit$fitted.values) %>%
      mutate(wt = if_else(A == 1, num/denom, (1-num)/(1-denom))) #Stabilized IPTW
    
    gset <- gset %>% mutate(wt2 = case_when(wt > quantile(gset$wt, 0.95) ~ quantile(gset$wt, 0.95),
                                            wt < quantile(gset$wt, 0.05) ~ quantile(gset$wt, 0.05),
                                            T ~ wt)) # Truncated IPTW
    # hist(denom_fit$fitted.values)
    # summary(denom_fit$fitted.values)
    # hist(gset$wt)
    # hist(gset$wt2)
    # plot(gset$A, gset$denom)
    # # mis-classification rate
    # sum(gset$A*(gset$denom<0.9999) + (1-gset$A)*(gset$denom>0.0001))/length(gset$A)
    # boxplot(gset$A*num/gset$denom)
    # summary(gset$A*num/gset$denom)
    # cbind(gset$A, gset$denom)
    # boxplot((1-gset$A)*num/(1-gset$denom))
    # summary((1-gset$A)*num/(1-gset$denom))
    # model_glm <- glm(data = gset, weight = wt2,
    #                  formula = as.formula(outForm),
    #                  family = "gaussian")
    require("survey")
    model_glm <-(svyglm(Y ~ A, design = svydesign(~ 1, weights = ~ wt2,
                                                  data = gset)))
    ATE <- summary(model_glm)$coefficients[2,1]
    library(sandwich)
    SE<-sqrt(diag(vcovHC(model_glm, type="HC0")))[2]
    # SE <- summary(model_glm)$coefficients[2,2]
    TYPE <- "GLM + IPW"
    ptm1 <- proc.time()
    t <- round((ptm1 - ptm0)[3],3)
    res <- rbind(res, cbind(ATE, SE, TYPE,t, no_cores))
  }
  
  if (doLASSO == 1){
    ptm0 <- proc.time()
    library(glmnet)
    # gset <- set1
    num <- summary(glm(data = gset, formula = A ~ 1, family = "gaussian"))$coefficient[1,1]
    x_exp <- model.matrix(as.formula(expForm), data=gset)[,-1]
    y_exp <- gset$A
    
    denom_fit <- cv.glmnet(y=y_exp, x=x_exp,
                           family = "binomial")
    # hist(expit(predict(denom_fit, newx = x_exp, s="lambda.min")))
    # quantile(expit(predict(denom_fit, newx = x_exp, s="lambda.min")), c(0.05,0.95))
    # fitted.values(denom_fit,s=lambda_min)
    # hist(denom_fit$fitted.values)
    # quantile(denom_fit$fitted.values, c(0.05,0.95))
    gset <- gset %>% mutate(denom = expit(predict(denom_fit, newx = x_exp, s="lambda.min"))) %>%
      mutate(wt = if_else(A == 1, num/denom, (1-num)/(1-denom))) #Stabilized IPTW
    # hist(gset$wt)
    #  quantile(gset$wt, c(0.05,0.95))
    
    gset <- gset %>% mutate(wt2 = case_when(wt > quantile(gset$wt, 0.95) ~ quantile(gset$wt, 0.95),
                                            wt < quantile(gset$wt, 0.05) ~ quantile(gset$wt, 0.05),
                                            T ~ wt)) # Truncated IPTW
    # plot(gset$A, gset$denom)
    # hist(gset$wt)
    # hist(gset$wt2)
    
    x <- model.matrix(as.formula(outForm), data=gset)[,-1]
    y <- gset$Y
    
    model_fit <- cv.glmnet(y=y, x=x,family = "gaussian")
    y.pred <- predict(model_fit, newx = x, s="lambda.min")
    # hist(y.pred)
    
    lambda_min <- model_fit$lambda.min
    coefficients <- as.matrix(coef(model_fit,s=lambda_min))
    coefficients <- coefficients[2:length(coefficients),]
    # coefficients[1] <- 1
    active.index <- which(coefficients != 0)
    # active.coefficients <- coefficients[active.index]

    x_min <- x[,active.index]
    # x_min <- x[,-1]
    # calculate propensity score, using min_lambda
    new_data <- data.frame(cbind(y,x_min))
    model_glm <- glm(y~.,data=new_data, weights=gset$wt2,family=gaussian())
    # hist(model_glm$fitted.values)
    # predict(model_fit, newx = x, s="lambda.min")
    # summary(model_fit, s="lambda.min")
    ATE <- summary(model_glm)$coefficients[2,1]
    SE <- summary(model_glm)$coefficients[2,2]
    TYPE <- "LASSO"
    ptm1 <- proc.time()
    t <- round((ptm1 - ptm0)[3],3)
    res <- rbind(res, cbind(ATE, SE, TYPE,t, no_cores))
  }
  
  # CV-TMLE default specification
  # nodes <- list(W = c("C1","C2","C3","C4","C5"), A = "A", Y = "Y")
  # tmle3_autofit <- tmle3(tmle_ATE(treatment_level = 1, control_level = 0), data.table::copy(tset), nodes, learner_list)
  # ATE <- tmle3_autofit$summary$tmle_est
  # SE <- tmle3_autofit$summary$se
  
  # GComputation
  if (doGComp == 1){
    # gset <- set1
    ptm0 <- proc.time()
    model_glm <- glm(data = gset,
                     formula = as.formula(outForm),
                     family = "gaussian")
    ate.vector <- rep(0,length(coef(model_glm)))
    ate.vector[2] <- 1; # ate.vector[8:12] <- colMeans(gset[,5+idx.handpick])
    
    ATE <- coef(model_glm) %*% ate.vector
    SE <- sqrt(t(ate.vector) %*% vcov(model_glm) %*% ate.vector)
    TYPE <- "GComp"
    ptm1 <- proc.time()
    t <- round((ptm1 - ptm0)[3],3)
    res <- rbind(res, cbind(ATE, SE, TYPE,t, no_cores))
  }
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
    ptm0 <- proc.time()
    outcome <- as.character(as.formula(outForm)[[2]])
    exposure <- as.character(as.formula(expForm)[[2]])
    TMLE <- TMLE(data=tset, exposure=exposure, outcome=outcome,
                 covarsT=all.vars(as.formula(expForm))[-1],
                 covarsO=all.vars(as.formula(outForm))[-1],
                 learners=aipw_lib,
                 control=control,
                 outForm=outForm, expForm=expForm)
    
    ATE <- TMLE$rd
    SE <- sqrt(TMLE$vd)
    TYPE <- "ManuTMLE"
    ptm1 <- proc.time()
    t <- round((ptm1 - ptm0)[3],3)
    res <- rbind(res, cbind(ATE, SE, TYPE,t, no_cores))
  }
  
  if (doManuTMLE_Auto==1){
    # Manual TMLE
    ptm0 <- proc.time()
    outcome <- as.character(as.formula(outForm)[[2]])
    exposure <- as.character(as.formula(expForm)[[2]])
    TMLE <- TMLE_Auto(data=tset, exposure=exposure, outcome=outcome,
                 covarsT=all.vars(as.formula(expForm))[-1],
                 covarsO=all.vars(as.formula(outForm))[-1],
                 learners=aipw_lib,
                 control=control)
    
    ATE <- TMLE$rd
    SE <- sqrt(TMLE$vd)
    TYPE <- "ManuTMLE_Auto"
    ptm1 <- proc.time()
    t <- round((ptm1 - ptm0)[3],3)
    res <- rbind(res, cbind(ATE, SE, TYPE,t, no_cores))
  }
  
  if (doDCTMLE==1){
    # DC-TMLE
    ptm0 <- proc.time()
    outcome <- as.character(as.formula(outForm)[[2]])
    exposure <- as.character(as.formula(expForm)[[2]])
    DCTMLE <- DCTMLE_Multiple(data=tset, exposure=exposure, outcome=outcome,
                              covarsT=all.vars(as.formula(expForm))[-1],
                              covarsO=all.vars(as.formula(outForm))[-1],
                              learners=aipw_lib,
                              control=control, num_cf=num_cf,
                              parallel=parallel)
    
    
    ATE <- DCTMLE$rd
    SE <- sqrt(DCTMLE$mvd)
    TYPE <- "DC-TMLE"
    ptm1 <- proc.time()
    t <- round((ptm1 - ptm0)[3],3)
    res <- rbind(res, cbind(ATE, SE, TYPE,t, no_cores))
  }
  
  
  return(res)
}