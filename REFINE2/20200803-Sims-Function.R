##########################################
# Plasmode simulation
##########################################
### Outcome and exposure with high-dim covars (N = 331 covars)

# Generate formula for low dim and high-dim, Y5/A1 and Y/A
make.formula <- function(outVar = "Y", expVar = "A", ver = "FULL", 
                         sims.ver = "plas",vars=c("VAR_1"),p=331,interact=F,
                         interact.w.exp = F,data.ver = "FULL"){
  if (sims.ver == "plas"){
    if (ver == "FULL"){
      outForm <- paste0(outVar," ~ ")
      if ( interact.w.exp == F){
      outForm <- paste0(outForm, expVar, " + ")
      }
      expForm <- paste0(expVar," ~ ")
      if (interact==T){
        outForm <- paste0(outForm," ( ")
        expForm <- paste0(expForm," (")
      }
      if ( interact.w.exp == T){
        outForm <- paste0(outForm, expVar)
      }
      # for(i in 1:length(vars)){
      p <- length(vars)
      for(i in 1:5){
        if(i == 1){
          expForm <- paste0(expForm, vars[i])
          if ( interact.w.exp == T){
            outForm <- paste0(outForm, " + ",vars[i])
          }else{
            outForm <- paste0(outForm, vars[i])
          }
            
        }else{expForm <- paste0(expForm, " + ", vars[i])
        outForm <- paste0(outForm, " + ",vars[i])
        }
      }
      if (interact==T){
        outForm <- paste0(outForm," )^2")
        expForm <- paste0(expForm," )^2")
      }
      if (p > 5){
        for(i in 6:p){
          expForm <- paste0(expForm, " + ", vars[i])
          outForm <- paste0(outForm, " + ",vars[i])
        }
      }
    }else{
      confVar <- "L0.a + L0.b + L0.c + L0.d + L0.e + L0.f + L0.g + L0.h + L0.i + L0.j + L0.k"
      outForm <- paste0(outVar," ~ ",expVar," + ", confVar)
      expForm <- paste0(expVar," ~ ", confVar)
    }

    return(list(outForm=outForm,
                expForm=expForm))
    # return(list(outForm=outForm,
    #             expForm=expForm,
    #             plas=plas))
  }
  else{
    expForm <- "A ~ C1 + C2 + C3 + C4 + C5"
    outForm <- "Y ~ A + C1 + C2 + C3 + C4 + C5"
    return(list(outForm=outForm,
                expForm=expForm))
  }

}

# Generate dataset. 
# Options: ver: high-dim vs low-dim
#         size, use.subset: N=1178 vs N < 1178; 
make.set <- function(ver = "FULL", size = 200, plas,use.subset, p=331){
  if (size < 1178) {
    use.subset = T
  }
  if (ver == "FULL"){
    plas <- data.frame(plas)[,c(1:(2+p),ncol(plas))]
  }
  else{
    ## Outcome and exposure model with original, restricted covariate set (N = 16 covars)
    plas_org <- plas %>%
      dplyr::select(A1,
                    L0.a = VAR_12, L0.b = VAR_4, L0.c = VAR_2, L0.d = VAR_11, L0.e = VAR_217, L0.f = VAR_13, L0.g = VAR_1,
                    L0.h = VAR_49, L0.i = VAR_42, L0.j = VAR_41, L0.k = VAR_3,
                    L1.a = VAR_9, L1.b = VAR_27, L1.c = VAR_28, L1.d = VAR_29,
                    L1.e = VAR_31, L1.f = VAR_16, L1.g = VAR_8, Y5)
    plas_org <- plas_org %>% add_column(id = c(1:nrow(plas_org))) %>% mutate(id = paste0("ID",id)) %>% data.frame(.)
    plas <- data.frame(plas_org)
  }
  if (use.subset==T){
    idxs <- sample(1:nrow(plas), size, replace = FALSE)
    plas <- plas[idxs,]
  }
  return(plas)
}

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
  # Adding error 
  Y <- Y + rnorm(Nsamp, 0, 4)
  I <- i
  id <- paste0("ID",1:Nsamp)
  if (sims.ver=="5var"){
    res <- c(as_tibble(cbind(C1, C2, C3, C4, C5, A, Y)))
    # res <- c(as_tibble(cbind(C1, C2, C3, C4, C5, A, Y, I,id)))
  }else{
    res <- c(as_tibble(cbind(C1, C2, C3, C4, C5, A, Y)))
  }
  
  return(res)
}


general.sim <- function(sims.ver = "plas",path="./plas_data.dta", estimateWithMore=T, plas.seed=1111,randVar=F,
                        isHandPick=T,idx.handpick=c(1,2,5,18, 217, 2:40),size=1178,
                        data.ver ="FULL", use.subset = F, generateA = T,
                        interact.w.exp = F, interact=T,plas_sim_N=500,
                        outForm=".", expForm="."){
  if (sims.ver == "plas") {
    # out_path <- paste0(path,"/Data/")
    # plas_org <- haven::read_dta(paste0("./plas_data.dta"))
    # plas_org <- haven::read_dta(path)
    plas_org <- read.csv(path, header=TRUE, stringsAsFactors=FALSE)
    # plas<- plas_org %>% add_column(id = c(1:nrow(plas_org))) %>% mutate(id = paste0("ID",id)) %>% data.frame(.)
    plas <- data.frame(plas_org)
    # rm(.Random.seed, envir=.GlobalEnv)
    # if (estimateWithMore == T){
    #   p.set <- 50
    # }else{
    #   p.set <- 100
    # } # This p.set only affect var.idx which only affects formula
    # 
    # set.seed(plas.seed)
    # if (randVar == T){
    #   var.idx <- sample(1:p.set, p.set, replace=FALSE)
    # }else{
    #   var.idx <- 1:p.set
    # }
    # 
    # if (isHandPick == T){
    #   var.idx <- idx.handpick
    # }
    
    
    # vars <- names(plas_org[3:333])[var.idx] # Exclude treatment and outcome
    # plas <- make.set(ver=data.ver, size = size, plas = plas_org, use.subset=use.subset,p=p.set)
    # plas <- make.set(ver=data.ver, size = size, plas = plas_org, use.subset=use.subset)
    
    # plas.formula <- make.formula("Y5", "A1", ver=data.ver,vars=vars, p=p.sim, interact=interact)
    # plas.formula <- make.formula("Y5", "A1", ver=data.ver,vars=vars, interact=interact,
    #                              interact.w.exp=interact.w.exp )
    
    # (outForm <- plas.formula$outForm)
    # expForm <- plas.formula$expForm
    # outForm <- outForm
    # expForm <- expForm
    
    
    # Effect_Size <- 6.6 # simulated risk difference = large change (e.g. absolute units)
    source("20200226-PlasmodeCont_Revised.R")
    plas_sims <- PlasmodeContNew(formulaOut = as.formula(outForm), objectOut = NULL,
                                 formulaExp = as.formula(expForm), objectExp = NULL,
                                 data = plas, idVar = "id", 
                                 effectOR = Effect_Size, MMOut = 1, MMExp = 1, 
                                 nsim = plas_sim_N, size = nrow(plas),
                                 exposedPrev = NULL,generateA=generateA,
                                 sims.ver=sims.ver)
    return(list(plas_sims=plas_sims,
                plas=plas,
                vars=vars,
                RegEff = plas_sims$RegEff))
  # }else if (sims.ver == "5var.then.plas"){
  #   set.seed(427)
  #   Effect_Size <- 6.6 
  #   plas <- data.frame(draw_sims(1))
  #   plas.formula <- make.formula("Y5", "A1", ver=data.ver,sims.ver = sims.ver)
  #   outForm <- plas.formula$outForm
  #   expForm <- plas.formula$expForm
  #   source("20200226-PlasmodeCont_Revised.R")
  #   set.seed(plas.seed)
  #   plas_sims <- PlasmodeContNew(formulaOut = as.formula(outForm), objectOut = NULL,
  #                                formulaExp = as.formula(expForm), objectExp = NULL,
  #                                data = plas, idVar = "id", 
  #                                effectOR = Effect_Size, MMOut = 1, MMExp = 1, 
  #                                nsim = Nsets, size = nrow(plas),
  #                                exposedPrev = NULL)
  #   return(list(plas_sims=plas_sims,
  #               plas=plas,
  #               vars=vars))
  }else{
    ##################################
    # Generate basic sim dataset
    ##################################
    # set.seed(plas.seed)
    # cat(plas.seed)
    
    sim_boots <- NULL
    # Effect_Size <- 6.6
    
    registerDoRNG(plas.seed)
    sim_boots <- foreach(i = 1:Nsets) %dopar% {
      draw_sims(i)
    }
    return(sim_boots)
  }
}
