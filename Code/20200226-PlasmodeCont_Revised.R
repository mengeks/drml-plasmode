# DEBUG. Initiate Within function variables
# plas <- make.set(ver=data.ver, size = size, plas = plas_org, use.subset=F)
# plas.formula <- make.formula("Y5", "A1", ver=data.ver)
# outForm <- plas.formula$outForm
# expForm <- plas.formula$expForm


# formulaOut = as.formula(outForm)
# objectOut = NULL
# formulaExp = as.formula(expForm)
# objectExp = NULL
# data = plas
# idVar = "id"
# effectOR = Effect_Size
# MMOut = 1
# MMExp = 1
# nsim = Nsets
# size = nrow(plas)
# exposedPrev = NULL


# set.seed(plas.seed)

PlasmodeContNew <- function(formulaOut=NULL, objectOut=NULL,formulaExp=NULL,objectExp=NULL, data, idVar,
                              effectOR =1, MMOut=1,MMExp=1, nsim, size, eventRate=NULL, exposedPrev=NULL,
                              outVar=NULL)
{
  require(glm2)
  require(mgcv)
  
  ########################################################################
  #### CONTINUOUS OUTCOME - Data set + formulas to estimate Outcome AND exposure.
  ########################################################################
  if(is.null(formulaOut)==FALSE & is.null(formulaExp)==FALSE & is.null(objectOut)==TRUE& is.null(objectExp)==TRUE)
  {
    outcome<- all.vars(formulaOut)[1]  ## selects the outcome variable
    exposure<- all.vars(formulaOut)[2] ## selects the exposure variable
    
    x <- data[order(data[,exposure]),] # order according to exposure status, unexposed first
    n <- nrow(x)
    n1 <- sum(as.numeric(x[,exposure]))    # number of exposed in real data
    n0 <- n - n1
    size1 <- round(ifelse(is.null(exposedPrev), n1, size*exposedPrev))  # desired number of exposed
    size0 <- size - size1 # desired number unexposed 
    if(size1 > n1 | size0 > n0) stop("Number of requested exposed or unexposed exceeds observed number -- reduce size")
    
    ## Estimate logit model for probability of exposure
    modExp<- glm2(formulaExp, family = "binomial", x, control=glm.control(trace=F))
    ## Design matrix used for exposure logistic regression
    XEXP<- gam(formulaExp, x, family = "binomial", fit = FALSE)$X
    
    # Adjusting the exposure prevalence in base cohort
    if(is.null(exposedPrev)) exposedPrev <- mean(modExp$y)
    bnewExp<- c(coef(modExp)[1], MMExp*coef(modExp)[-1])
    
      XbnewExp<- as.vector(XEXP%*%bnewExp)
      fnExp<- function(d)mean(plogis(d+XbnewExp))-exposedPrev
      deltaExp <- uniroot(fnExp, lower=-20, upper=20)$root
      Probexp <- plogis(deltaExp+XbnewExp)
      # resample new exposure measures
      if (generateA==T){
        x[exposure] <- rbinom(size,1,Probexp)
      }
      rm(modExp, XEXP)
    
    # Compute outcome model, using new exposures
    modOutCont <- glm2(formulaOut, family = "gaussian", x, control=glm.control(trace=F)) ## outcome model coefficients
    X <- gam(formulaOut, x, family = "gaussian", fit = FALSE)$X ## extract outcome model matrix
    
    bnew <- c(coef(modOutCont)[1], MMOut*coef(modOutCont)[-1]) # find intercept value needed to get approximate mean outcome
    bnew <- replace(bnew, names(coef(modOutCont)) == exposure, effectOR) # "inject" desired RD (beta coeff)
    EYnew <- as.vector(X %*% bnew)
    # err <- rnorm(nrow(data),0,var(residuals(modOutCont)))
    # Xbnew <- EYnew + err# generate new outcome measure with noise
    
    # 7 Aug 2020: We add residuals randomly after 
    err <- residuals(modOutCont)
    Xbnew <- EYnew
    
    # # draw bootstrap sets
    ids <- ynew <- expnew <- data.frame(matrix(nrow = size, ncol = nsim)) # intialize replacement matrices 
    p1.vec <- p0.vec <- RR <- RD <- vector('numeric', length = nsim)
    As <- vector('numeric', length = nsim)#DEBUG
    # sim<-1 #DEBUG
    for(sim in 1:nsim) {
      if(sim == 1 | (sim %% (nsim/2)) == 0) print(paste0("Drawing sample: ", sim, " (of ", nsim, ") ..."))
      idxs <- sample(1:n, size, replace = TRUE) # resample from all rows, with replacement
      if (sims.ver=="plas"){
        ids[1:size,sim] <- x[idxs, idVar] # specify the exact ID for those rows
      }else{
        ids[1:size,sim] <-"ID1" # DUMMY
      }
      expnew[,sim]<- x[idxs, exposure] # draw exposure values 
      # ynew[,sim] <- Xbnew[idxs] + err # No bootstrap on errors
      ynew[,sim] <- Xbnew[idxs] + err[sample(1:n, size, replace = TRUE)] # Bootstrap on errors
      # ynew[,sim] <- Xbnew[idxs] # no error terms
      datasim<-X[idxs,]
      datasim[,2]<-1
      p_1<- as.vector(datasim %*% bnew)
      datasim[,2]<-0
      p_0<- as.vector(datasim %*% bnew)
      RR[sim]<-mean(p_1)/mean(p_0)
      RD[sim]<-mean(p_1)-mean(p_0)
      p1.vec[sim]<-mean(p_1)
      p0.vec[sim]<-mean(p_0)
    }
    ARR<-mean(RR); ARD<-mean(RD)
    Ap1 <- mean(p1.vec); Ap0 <-mean(p0.vec)
    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(ynew) <- paste("OUTCOME", 1:nsim, sep = "")
    names(expnew)<-paste("EXPOSURE",1:nsim, sep = "")
    sim_out_bin<-data.frame(ids, ynew,expnew)
    
    return(list(TrueOutBeta = bnew, TrueExpBeta = bnewExp, 
                RR=ARR,RD=ARD,Sim_Data = sim_out_bin,
                outForm=formulaOut, expForm=formulaExp,
                p1=Ap1, p0=Ap0))
  }
  
  ########################################################################
  #### CONTINUOUS OUTCOME - Data set + formulas to estimate Exposure ONLY. -- outVar must be specified
  ########################################################################
  else if(is.null(formulaOut)==TRUE & is.null(formulaExp)==FALSE & is.null(objectOut)==TRUE& is.null(objectExp)==TRUE)
  {
    
    if(is.null(outVar)) stop("Outcome variable (outVar) must be specified if only exposure models are provided.")
    
    exposure <- all.vars(formulaExp)[1]
    #exposure <- as.formula(formulaExp)[[2]] ##selects the exposure variable
    
    x <- data[order(data[,exposure]),] # order according to exposure status, unexposed first
    n <- nrow(x)
    n1 <- sum(x[,exposure])     # number of exposed in real data
    n0 <- n - n1
    size1 <- round(ifelse(is.null(exposedPrev), n1, size*exposedPrev))  # desired number of exposed
    size0 <- size - size1
    if(size1 > n1 | size0 > n0) stop("Number of requested exposed or unexposed exceeds observed number -- reduce size")
    
    ## Estimate logit model for probability of exposure
    modExp<- glm2(formulaExp, family = "binomial", x, control=glm.control(trace=F))
    ## Design matrix used for exposure logistic regression
    XEXP<- gam(formulaExp, x, family = "binomial", fit = FALSE)$X ### corrected again here - JYH
    # Adjusting the exposure prevalence in base cohort
    if(is.null(exposedPrev))exposedPrev<- mean(modExp$y)
    bnewExp<- c(coef(modExp)[1], MMExp*coef(modExp)[-1])
    XbnewExp<- as.vector(XEXP%*%bnewExp)
    fnExp<- function(d)mean(plogis(d+XbnewExp))-exposedPrev
    deltaExp<- uniroot(fnExp, lower=-20, upper=20)$root
    Probexp<- plogis(deltaExp+XbnewExp)
    # resample new exposure measures
    x[exposure] <- rbinom(size,1,Probexp)  
    rm(modExp, XEXP)
    
    # compute outcomes using new exposure
    modOutCont <- glm2(as.formula(paste0(outVar, " ~ ", exposure)), family = "gaussian", data = x, control=glm.control(trace=F)) ## UNWEIGHTED - outcome model coefficients
    X <- gam(as.formula(paste0(outVar, " ~ ", exposure)), family = "gaussian", data = x, fit = FALSE)$X ## UNWEIGHTED - extract outcome model matrix
    
    bnew <- c(coef(modOutCont)[1], MMOut*coef(modOutCont)[-1]) # find intercept value needed to get approximate mean outcome
    bnew <- replace(bnew, names(coef(modOutCont)) == exposure, effectOR) # "inject" desired RD (beta coeff)
    Xbnew <- as.vector(X %*% bnew)+rnorm(nrow(data),0,var(residuals(modOutCont))) # generate new outcome measure with noise
  
    ids <- ynew <- expnew<-data.frame(matrix(nrow = size, ncol = nsim))
    for(sim in 1:nsim) {
      if(sim == 1 | (sim %% (nsim/2)) == 0) print(paste0("Drawing sample: ", sim, " (of ", nsim, ") ..."))
      idxs <- sample(1:n, size, replace = TRUE) # sample unexposed (located in rows 1:n0 of x)
      ids[1:size,sim] <- x[idxs, idVar]
      expnew[,sim]<- x[idxs, exposure]
      ynew[,sim] <- Xbnew[idxs] 
    }
    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(ynew) <- paste("OUTCOME", 1:nsim, sep = "")
    names(expnew)<-paste("EXPOSURE",1:nsim, sep = "")
    sim_out_bin<-data.frame(ids,ynew,expnew)
    
    return(list(TrueExpBeta = bnewExp, TrueOutBeta = bnew, Sim_Data = sim_out_bin))
  }
  
  ########################################################################
  #### CONTINUOUS OUTCOME - Data set + formulas to estimate Outcome ONLY.
  ########################################################################
  else if(is.null(formulaOut)==FALSE & is.null(formulaExp)==TRUE & is.null(objectOut)==TRUE& is.null(objectExp)==TRUE)
  {
    outcome<- all.vars(formulaOut)[1] ## selects the outcome variable
    exposure<- all.vars(formulaOut)[2] ##selects the exposure variable
    
    x <- data[order(data[,exposure]),] # order according to exposure status, unexposed first
    n <- nrow(x)
    n1 <- sum(x[,exposure])     # number of exposed in real data
    n0 <- n - n1
    size1 <- round(ifelse(is.null(exposedPrev), n1, size*exposedPrev))  # desired number of exposed
    size0 <- size - size1
    if(size1 > n1 | size0 > n0) stop("Number of requested exposed or unexposed exceeds observed number -- reduce size")
    
    # Set / adjust exposure prevalence and resample
    if(is.null(exposedPrev))exposedPrev<- sum(x[exposure])/nrow(x)
    x[exposure] <- rbinom(size,1,exposedPrev)  
    
    # compute outcomes using outcome model
    modOutCont <- glm2(formulaOut, family = "gaussian", data = x, control=glm.control(trace=F)) ## UNWEIGHTED - outcome model coefficients
    X <- gam(formulaOut, family = "gaussian", data = x, fit = FALSE)$X ## UNWEIGHTED - extract outcome model matrix
    
    bnew <- c(coef(modOutCont)[1], MMOut*coef(modOutCont)[-1]) # find intercept value needed to get approximate mean outcome
    bnew <- replace(bnew, names(coef(modOutCont)) == exposure, effectOR) # "inject" desired RD (beta coeff)
    Xbnew <- as.vector(X %*% bnew)+rnorm(nrow(data),0,var(residuals(modOutCont))) # generate new outcome measure with noise
    
    # draw bootstrap sets
    ids <- ynew <- expnew <- data.frame(matrix(nrow = size, ncol = nsim)) # intialize replacement matrices 
    RR <- RD <- vector('numeric', length = nsim)
    for(sim in 1:nsim) {
      if(sim == 1 | (sim %% (nsim/2)) == 0) print(paste0("Drawing sample: ", sim, " (of ", nsim, ") ..."))
      idxs <- sample(1:n, size, replace = TRUE) # resample from all rows, with replacement
      ids[1:size,sim] <- x[idxs, idVar] # specify the exact ID for those rows
      expnew[,sim]<- x[idxs, exposure]
      ynew[,sim] <- Xbnew[idxs] 
      datasim<-X[idxs,]
      datasim[,2]<-1
      p_1<- as.vector(datasim %*% bnew)
      datasim[,2]<-0
      p_0<- as.vector(datasim %*% bnew)
      RR[sim]<-mean(p_1)/mean(p_0)
      RD[sim]<-mean(p_1)-mean(p_0)
    }
    ARR<-mean(RR)
    ARD<-mean(RD)
    
    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(ynew) <- paste("OUTCOME", 1:nsim, sep = "")
    names(expnew)<-paste("EXPOSURE",1:nsim, sep = "")
    sim_out_bin<-data.frame(ids,ynew,expnew)
  
  return(list(TrueOutBeta = bnew, RR = ARR, RD = ARD, Sim_Data = sim_out_bin))
  }
  
  ########################################################################
  #### CONTINUOUS OUTCOME - Regression objects for Exposure AND Outcome.
  ########################################################################
  else if(is.null(formulaOut)==TRUE & is.null(formulaExp)==TRUE & is.null(objectOut)==FALSE & is.null(objectExp)==FALSE)
  {
    DesMatOut<- model.matrix(objectOut)
    exposure<- all.vars(objectOut$formula)[2] # selects the exposure variable
    outcome <- all.vars(objectOut$formula)[1] # selects the outcome variable
    dataOut<- cbind(objectOut$data[idVar], 
                    objectOut$data[outcome], 
                    as.data.frame(DesMatOut)) ## use as basis for final estimations for simplicity -- JYH
  
    n <- nrow(dataOut)
    n1 <- sum(dataOut[,exposure])     # number of exposed in real data
    n0 <- n - n1
    size1 <- round(ifelse(is.null(exposedPrev), n1, size*exposedPrev))  # desired number of exposed
    size0 <- size - size1
    if(size1 > n1 | size0 > n0) stop("Number of requested exposed or unexposed exceeds observed number -- reduce size")
    
    # Set / adjust exposure prevalence
    if(is.null(exposedPrev))exposedPrev<- sum(dataOut[exposure])/nrow(dataOut)
    DesMatExp<-model.matrix(objectExp) ### Assume rows will have to be sorted identically to outcome object - JYH
    XEXP<- as.data.frame(DesMatExp)
    ExpCoeff<-coef(objectExp)
    bnewExp<- c(ExpCoeff[1], MMExp*ExpCoeff[-1])
    XbnewExp<- as.vector(bnewExp%*%t(XEXP))
    fnExp<- function(d)mean(plogis(d+XbnewExp))-exposedPrev
    deltaExp<- uniroot(fnExp, lower=-20, upper=20)$root
    Probexp<- plogis(deltaExp+XbnewExp)
    
    # resample new exposure measures
    dataOut[exposure] <- rbinom(size,1,Probexp) 
    
    # replace exposure values with newly sampled exposure and re-estimate outcome model
    modOutCont <- glm2::glm2(objectOut$formula, family = "gaussian", dataOut, control=glm.control(trace=F)) ## outcome model coefficients
    X <- mgcv::gam(objectOut$formula, dataOut, family = "gaussian", fit = FALSE)$X ## extract outcome model matrix
    
    # Set new intercept value to get approximate desired event rate under new parameters
    bnew <- c(coef(modOutCont)[1], MMOut*coef(modOutCont)[-1]) # find intercept value needed to get approximate mean outcome
    bnew <- replace(bnew, names(coef(modOutCont)) == exposure, effectOR) # "inject" desired RD (beta coeff)
    Xbnew <- as.vector(X %*% bnew)+rnorm(nrow(dataOut),0,var(residuals(modOutCont))) # generate new outcome measure with noise

    # draw bootstrap sets
    ids <- ynew <- expnew <- data.frame(matrix(nrow = size, ncol = nsim)) # intialize replacement matrices 
    RR <- RD <- vector('numeric', length = nsim)
    for(sim in 1:nsim) {
      if(sim == 1 | (sim %% (nsim/2)) == 0) print(paste0("Drawing sample: ", sim, " (of ", nsim, ") ..."))
      idxs <- sample(1:n, size, replace = TRUE) # resample from all rows, with replacement
      ids[1:size,sim] <- dataOut[idxs, idVar] # specify the exact ID for those rows
      expnew[,sim]<- dataOut[idxs, exposure] # draw exposure values 
      ynew[,sim] <- Xbnew[idxs] ### draw respective outcome values 
      datasim<-X[idxs,]
      datasim[,2]<-1
      p_1<- as.vector(datasim %*% bnew)
      datasim[,2]<-0
      p_0<- as.vector(datasim %*% bnew)
      RR[sim]<-mean(p_1)/mean(p_0)
      RD[sim]<-mean(p_1)-mean(p_0)
    }
    ARR<-mean(RR)
    ARD<-mean(RD)
    
    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(ynew) <- paste("OUTCOME", 1:nsim, sep = "")
    names(expnew)<-paste("EXPOSURE",1:nsim, sep = "")
    sim_out_bin<-data.frame(ids, ynew,expnew)
    
    return(list(TrueOutBeta = bnew, TrueExpBeta = bnewExp, RR=ARR,RD=ARD,Sim_Data = sim_out_bin))
  }
 
  ########################################################################
  #### CONTINUOUS OUTCOME - Regression objects for Exposure ONLY. -- outVar must be specified
  ########################################################################
  else if(is.null(formulaOut)==TRUE & is.null(formulaExp)==TRUE & is.null(objectOut)==TRUE & is.null(objectExp)==FALSE)
  {
    if(is.null(outVar)) stop("Outcome variable (outVar) must be specified if only exposure models are provided.")
    
    exposure <- all.vars(objectExp$formula)[1] # select exposure variable
    DesMatExp<-model.matrix(objectExp) 
    dataOut<- cbind(objectExp$data[idVar],
                    objectExp$data[exposure],
                    objectExp$data[outVar],  ## Thus, exposure model object needs to be fit in the dataset with the outcome variable
                    as.data.frame(DesMatExp)) ## use as basis for final estimations for simplicity -- JYH
    
    n <- nrow(dataOut)
    n1 <- sum(dataOut[,exposure])     # number of exposed in real data
    n0 <- n - n1
    size1 <- round(ifelse(is.null(exposedPrev), n1, size*exposedPrev))  # desired number of exposed
    size0 <- size - size1
    if(size1 > n1 | size0 > n0) stop("Number of requested exposed or unexposed exceeds observed number -- reduce size")
    
    # Set / adjust exposure prevalence
    if(is.null(exposedPrev))exposedPrev<- sum(dataOut[exposure])/nrow(dataOut)
    XEXP<- as.data.frame(DesMatExp)
    ExpCoeff<-coef(objectExp)
    bnewExp<- c(ExpCoeff[1], MMExp*ExpCoeff[-1])
    XbnewExp<- as.vector(bnewExp%*%t(XEXP))
    fnExp<- function(d)mean(plogis(d+XbnewExp))-exposedPrev
    deltaExp<- uniroot(fnExp, lower=-20, upper=20)$root
    Probexp<- plogis(deltaExp+XbnewExp)
    
    # resample new exposure measures
    dataOut[exposure] <- rbinom(size,1,Probexp) 
    
    # replace exposure values with newly sampled exposure and re-estimate outcome model
    modOutCont <- glm2::glm2(as.formula(paste0(outVar, " ~ ", exposure)), family = "gaussian", dataOut, control=glm.control(trace=F)) ## outcome model coefficients
    X <- mgcv::gam(as.formula(paste0(outVar, " ~ ", exposure)), dataOut, family = "gaussian", fit = FALSE)$X ## extract outcome model matrix
    
    # Set new intercept value to get approximate desired event rate under new parameters
    bnew <- c(coef(modOutCont)[1], MMOut*coef(modOutCont)[-1]) # find intercept value needed to get approximate mean outcome
    bnew <- replace(bnew, names(coef(modOutCont)) == exposure, effectOR) # "inject" desired RD (beta coeff)
    Xbnew <- as.vector(X %*% bnew)+rnorm(nrow(dataOut),0,var(residuals(modOutCont))) # generate new outcome measure with noise
    
    # draw bootstrap sets
    ids <- ynew <- expnew <- data.frame(matrix(nrow = size, ncol = nsim)) # intialize replacement matrices 
    RR <- RD <- vector('numeric', length = nsim)
    for(sim in 1:nsim) {
      if(sim == 1 | (sim %% (nsim/2)) == 0) print(paste0("Drawing sample: ", sim, " (of ", nsim, ") ..."))
      idxs <- sample(1:n, size, replace = TRUE) # resample from all rows, with replacement
      ids[1:size,sim] <- dataOut[idxs, idVar] # specify the exact ID for those rows
      expnew[,sim]<- dataOut[idxs, exposure] # draw exposure values 
      ynew[,sim] <- Xbnew[idxs] ### draw respective outcome values 
      datasim<-X[idxs,]
      datasim[,2]<-1
      p_1<- as.vector(datasim %*% bnew)
      datasim[,2]<-0
      p_0<- as.vector(datasim %*% bnew)
      RR[sim]<-mean(p_1)/mean(p_0)
      RD[sim]<-mean(p_1)-mean(p_0)
    }
    ARR<-mean(RR)
    ARD<-mean(RD)
    
    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(ynew) <- paste("OUTCOME", 1:nsim, sep = "")
    names(expnew)<-paste("EXPOSURE",1:nsim, sep = "")
    sim_out_bin<-data.frame(ids, ynew,expnew)
    
    return(list(TrueOutBeta = bnew, TrueExpBeta = bnewExp, RR=ARR,RD=ARD,Sim_Data = sim_out_bin))
  } 
  
  ########################################################################
  #### CONTINUOUS OUTCOME - Regression objects for Outcome ONLY.
  ########################################################################
  else if(is.null(formulaOut)==TRUE & is.null(formulaExp)==TRUE & is.null(objectOut)==FALSE & is.null(objectExp)==TRUE)
  {
    DesMatOut<- model.matrix(objectOut)
    exposure<- all.vars(objectOut$formula)[2] # selects the exposure variable
    outcome <- all.vars(objectOut$formula)[1] # selects the outcome variable
    dataOut<- cbind(objectOut$data[idVar], 
                    objectOut$data[outcome], 
                    as.data.frame(DesMatOut)) ## use as basis for final estimations for simplicity -- JYH
    
    n <- nrow(dataOut)
    n1 <- sum(dataOut[,exposure])     # number of exposed in real data
    n0 <- n - n1
    size1 <- round(ifelse(is.null(exposedPrev), n1, size*exposedPrev))  # desired number of exposed
    size0 <- size - size1
    if(size1 > n1 | size0 > n0) stop("Number of requested exposed or unexposed exceeds observed number -- reduce size")
    
    # Set / adjust exposure prevalence and resample
    if(is.null(exposedPrev))exposedPrev<- sum(dataOut[exposure])/nrow(dataOut)
    dataOut[exposure] <- rbinom(size,1,exposedPrev) 
    
    # replace exposure values with newly sampled exposure and re-estimate outcome model
    modOutCont <- glm2::glm2(objectOut$formula, family = "gaussian", dataOut, control=glm.control(trace=F)) ## outcome model coefficients
    X <- mgcv::gam(objectOut$formula, dataOut, family = "gaussian", fit = FALSE)$X ## extract outcome model matrix
    
    # Set new intercept value to get approximate desired event rate under new parameters
    bnew <- c(coef(modOutCont)[1], MMOut*coef(modOutCont)[-1]) # find intercept value needed to get approximate mean outcome
    bnew <- replace(bnew, names(coef(modOutCont)) == exposure, effectOR) # "inject" desired RD (beta coeff)
    Xbnew <- as.vector(X %*% bnew)+rnorm(nrow(dataOut),0,var(residuals(modOutCont))) # generate new outcome measure with noise
    
    # draw bootstrap sets
    ids <- ynew <- expnew <- data.frame(matrix(nrow = size, ncol = nsim)) # intialize replacement matrices 
    RR <- RD <- vector('numeric', length = nsim)
    for(sim in 1:nsim) {
      if(sim == 1 | (sim %% (nsim/2)) == 0) print(paste0("Drawing sample: ", sim, " (of ", nsim, ") ..."))
      idxs <- sample(1:n, size, replace = TRUE) # resample from all rows, with replacement
      ids[1:size,sim] <- dataOut[idxs, idVar] # specify the exact ID for those rows
      expnew[,sim]<- dataOut[idxs, exposure] # draw exposure values 
      ynew[,sim] <- Xbnew[idxs] ### draw respective outcome values 
      datasim<-X[idxs,]
      datasim[,2]<-1
      p_1<- as.vector(datasim %*% bnew)
      datasim[,2]<-0
      p_0<- as.vector(datasim %*% bnew)
      RR[sim]<-mean(p_1)/mean(p_0)
      RD[sim]<-mean(p_1)-mean(p_0)
    }
    ARR<-mean(RR)
    ARD<-mean(RD)
    
    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(ynew) <- paste("OUTCOME", 1:nsim, sep = "")
    names(expnew)<-paste("EXPOSURE",1:nsim, sep = "")
    sim_out_bin<-data.frame(ids, ynew,expnew)
    
    return(list(TrueOutBeta = bnew, RR=ARR,RD=ARD,Sim_Data = sim_out_bin))
  }
}



#########################################
#########################################
#########################################
# ## TEST BY LOOPING THROUGH SIMULATIONS
# library(tidyverse)
# plas <- readxl::read_xlsx("sim_data.xlsx") %>% as.data.frame(.)
# mean(plas[which(plas$A1==1),]$Y5) - mean(plas[which(plas$A1==0),]$Y5)
# #plas <- readxl::read_excel("sim_data.xlsx")
# 
# objectOut <- glm(formula = Y5 ~ A1 + L0.a + L0.d + L0.e + L1.a, family = gaussian, data = plas)
# objectExp <- glm(formula = A1 ~ L0.a + L0.d + L0.e + L1.a, family = binomial, data = plas)
# formulaOut <- as.formula("Y5 ~ A1 + L0.a + L0.d + L0.e + L1.a")
# formulaExp <- as.formula("A1 ~ L0.a + L0.d + L0.e + L1.a")
# 
# {
#   # Generate sims for all scenarios
#   Effect_Size = 6.6
#   N_sims = 500
#   set.seed(12345)
#   print("form_out_exp")
#   form_out_exp <- PlasmodeContNew(formulaOut = formulaOut, formulaExp = formulaExp,
#                               data = plas, idVar = "id",
#                               MMExp = 1, MMOut = 1,
#                               effectOR = Effect_Size, nsim = N_sims,
#                              size = nrow(plas))
#   print("form_out")
#   form_out <- PlasmodeContNew(formulaOut = formulaOut, formulaExp = NULL,
#                            data = plas, idVar = "id",
#                            MMExp = 1, MMOut = 1,
#                            effectOR = Effect_Size, nsim = N_sims,
#                            size = nrow(plas))
#   print("form_exp")
#   form_exp <- PlasmodeContNew(formulaOut = NULL, formulaExp = formulaExp,
#                               data = plas, idVar = "id",
#                               MMExp = 1, MMOut = 1,
#                               effectOR = Effect_Size, nsim = N_sims,
#                               size = nrow(plas), outVar = "Y5")
#   print("obj_out_exp")
#   obj_out_exp <- PlasmodeContNew(objectOut = objectOut, objectExp = objectExp,
#                               data = NULL, idVar = "id",
#                               MMExp = 1, MMOut = 1,
#                               effectOR = Effect_Size, nsim = N_sims,
#                               size = nrow(plas))
#   print("obj_out")
#   obj_out <- PlasmodeContNew(objectOut = objectOut, objectExp = NULL,
#                                  data = NULL, idVar = "id",
#                                  MMExp = 1, MMOut = 1,
#                                  effectOR = Effect_Size, nsim = N_sims,
#                                  size = nrow(plas))
#   print("obj_exp")
#   obj_exp <- PlasmodeContNew(objectOut = NULL, objectExp = objectExp,
#                              data = NULL, idVar = "id",
#                              MMExp = 1, MMOut = 1,
#                              effectOR = Effect_Size, nsim = N_sims,
#                              size = nrow(plas), outVar = "Y5")
# }
# 
# 
# # estimate ate, se, bias, mse, coverage for each type
# sim_names <- c("form_out_exp", "form_out", "form_exp", "obj_out_exp", "obj_out", "obj_exp")
# sim_res <- NULL
# 
# 
# for(t in 1:length(sim_names)){
#   print(sim_names[t])
#   sdf <- NULL
#   sims1 <- eval(as.name(sim_names[t]))
#   for(i in 1:N_sims){
#     stmp <- cbind(id = sims1$Sim_Data[i],
#                   A1 = sims1$Sim_Data[i + 2*N_sims],
#                   Y5 = sims1$Sim_Data[i + N_sims])
#     colnames(stmp) <- c("id", "A1", "Y5")
#     suppressMessages(stmp <- left_join(as_tibble(stmp), dplyr::select(plas, -Y5, -A1)))
#     # mean(stmp[which(stmp$A1==1),]$Y5) - mean(stmp[which(stmp$A1==0),]$Y5)
#     if(i%%50 == 0){print(paste0("Estimating set: ", i, " (out of ", N_sims, ")"))}
#     pn <- summary(glm(data = stmp, formula = A1 ~ 1, family = "gaussian"))$coefficient[1,1]
#     if(sim_names[t] %in% c("form_out_exp", "form_exp")){pd <- glm(data = stmp, formula = formulaExp, family = "binomial")}
#     if(sim_names[t] %in% c("obj_out_exp", "obj_exp")){pd <- glm(data = stmp, formula = objectExp$formula, family = "binomial")}
#     stmp <- stmp %>% add_column(pd = pd$fitted.values) %>%
#       mutate(wt0 = 1, # unweighted
#              wt1 = if_else(A1 == 1, pn/pd, (1-pn)/(1-pd))) %>% # Stabilized IPTW
#       mutate(wt2 = case_when(wt1 > quantile(wt1, 0.95) ~ quantile(wt1, 0.95),
#                              wt1 < quantile(wt1, 0.05) ~ quantile(wt1, 0.05),
#                              T ~ wt1)) # Truncated IPTW
#     if(sim_names[t] %in% c("obj_exp", "form_exp")){sres <- stmp %>% glm(formula = Y5 ~ A1, family = "gaussian", weight = wt2) %>% summary()}
#     if(sim_names[t] %in% c("obj_out_exp")){sres <- stmp %>% glm(formula = objectOut$formula, family = "gaussian", weight = wt2) %>% summary()}
#     if(sim_names[t] %in% c("obj_out")){sres <- stmp %>% glm(formula = objectOut$formula, family = "gaussian", weight = wt0) %>% summary()}
#     if(sim_names[t] %in% c("form_out_exp")){sres <- stmp %>% glm(formula = formulaOut, family = "gaussian", weight = wt2) %>% summary()}
#     if(sim_names[t] %in% c("form_out")){sres <- stmp %>% glm(formula = formulaOut, family = "gaussian", weight = wt0) %>% summary()}
#     BETA <- sres$coefficients[2,1]
#     SE <- sres$coefficients[2,2]
#     I <- i
#     sdf <- rbind(sdf, cbind(BETA, SE, I))
#   }
#   sdf <- as_tibble(sdf) %>% mutate(BETA = BETA, SE = SE,
#                                    LL = BETA - 1.96*SE, UL = BETA + 1.96*SE,
#                                    mu = mean(BETA), bias = as.double(BETA) - Effect_Size)
# 
#   # pt <- sdf %>% ggplot() + geom_errorbar(aes(x = I, ymin = LL, ymax = UL)) +
#   #   geom_hline(aes(yintercept = mu), linetype = "dashed") +
#   #   geom_hline(aes(yintercept = Effect_Size), linetype = "solid")
#   # plot(pt)
#   # sdf %>% summarize(mu_ATE = mean(BETA), med_ATE = median(BETA),
#   #                   mu_SE = mean(SE), med_SE = median(SE),
#   #                   mu_bias = mean(bias), med_bias = median(bias),
#   #                   var = var(BETA), MSE = var + mu_bias^2,
#   #                   coverage = sum((LL <= Effect_Size & UL >= Effect_Size)/N_sims))
# 
#   res <- sdf %>% summarize(mu_ATE = mean(BETA), mu_SE = mean(SE),
#                     mu_bias = mean(bias), MSE = var(BETA) + mu_bias^2,
#                     coverage = sum((LL <= Effect_Size & UL >= Effect_Size)/N_sims)) %>% as.numeric()
#   sim_res <- rbind(sim_res, res)
# }
# colnames(sim_res) <- c("BETA", "SE", "BIAS", "MSE", "COVERAGE")
# as_tibble(sim_res) %>% add_column(sim_names)
