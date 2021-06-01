##################################
## SUMMARIZE AND VISUALIZE
##################################
library(ggplot2)
library(tidyverse)
  setwd("~/Desktop/HuangGroup/cvtmle_plasmode/Code")
  # Effect_Size <- 6.6
  # Effect_Size <- 5.51788
  Effect_Size <- 1.190723
  
  path <- "~/Desktop/HuangGroup/cvtmle_plasmode/Code/cluster/no2/"
  folder <- "RDataFiles";
  # folder<-"ResBefore5Apr2021"
  is.timing <- T
  tim.str <- ifelse(is.timing,"-timing","")
  # mtd.lst <- c("IPW","GComp","TMLE","AIPW","DCTMLE","DCAIPW")
  mtd.lst  <- c("IPW","TMLE","AIPW","DCTMLE","DCAIPW")
  # mtd.lst <- c("GComp","TMLE","AIPW","DCTMLE","DCAIPW")
  # mtd.lst  <- c("TMLE","AIPW","DCTMLE","DCAIPW")
  
  adjustPath <- function(res.path){
    adjusted=F
    tryCatch(
      {
        message(paste0("Try loading ", res.path))
        load(res.path)
      },
      error=function(cond) {
        message(paste("There is no file at", res.path))
        res.path <- paste0(path,folder, "/result-",mtd,"-par",".RData")
        adjusted=T
        message(paste("Set path to the non-timing version:", res.path))
        return(NULL)
      },
      finally={
        a <- 1
      }
    )
    return(list(res.path=res.path, adjusted=adjusted))
  }
  
  is.timing <- T;tim.str <- ifelse(is.timing,"-timing","")
  {
    out.lst <- c("med_bias","med_SE","coverage","var", "med_t")
    out.str.lst <- character(length(out.lst))
    source("20200816-Result-Summary.R")
    i <- 1
    # load(paste0(path,folder, "/result-GComp","-par",tim.str,".RData"))
    for (out.var in out.lst){
      out.table.str <- ""
      for (mtd in mtd.lst){
        # if ((mtd  != "IPW") || (mtd  != "GComp")){
        #   is.timing <- T;tim.str <- ifelse(is.timing,"-timing","")
        #   res.path <- paste0(path,folder, "/result-",mtd,"-non-par",tim.str,".RData")
        # }else{
          res.path <- paste0(path,folder, "/result-",mtd,"-par",tim.str,".RData")
          load(res.path)
          path.obj <- adjustPath(res.path)
          res.path <- path.obj$res.path
          out.table.str <- append.res(out.table.str, res.path, out.var, path.obj$adjusted)

        if ((mtd  != "IPW") && (mtd  != "GComp")){ # append non-par result
          res.path <- paste0(path,folder, "/result-",mtd,"-non-par",tim.str,".RData")
          path.obj <- adjustPath(res.path)
          res.path <- path.obj$res.path
          out.table.str <- append.res(out.table.str, res.path, out.var, path.obj$adjusted)
        }
      }
      out.str.lst[i] <- out.table.str
      i <- i+1
    }
    out.str.lst[1] <- paste0("Bias ($\\times$ 100)",out.str.lst[1], "\\","\\  [2pt] \n")
    out.str.lst[2] <- paste0("SE ",out.str.lst[2], "\\","\\  [2pt] \n")
    out.str.lst[3] <- paste0("CI covg. ",out.str.lst[3], "\\","\\  [2pt] \n")
    out.str.lst[4] <- paste0("BVar ",out.str.lst[4], "\\","\\  [10pt] \n")
    out.str.lst[5] <- paste0("30 cores ",out.str.lst[5], "\\","\\  [2pt] \n")
  }

  cat(out.str.lst)
  
  
 






  {
    {
      source("20200816-Result-Summary.R")
      path.lst <- c("~/Desktop/HuangGroup/cvtmle_plasmode/Code/cluster/no4/")
      situ.lst <- c("A.bad")
      Effect_Size = 6.6
      mtd.lst <- c("IPW","TMLE","AIPW","DCTMLE","DCAIPW")
      factor.lst <- c("IPW","TMLE.par","TMLE.non.par","AIPW.par","AIPW.non.par",
                      "DC-TMLE.par","DC-TMLE.non.par", "DC-AIPW.par","DC-AIPW.non.par")
      violin.and.bar.plot.by.situ(path.lst, situ.lst, ci.axis.coeff=0.2)
    }
    ggsave("/Users/garethalex/Desktop/HuangGroup/cvtmle_plasmode/Output/situA.png",
           width = 14, height = 14, dpi = 300, units = "in", device='png')
  }
  
  {
    {
      source("20200816-Result-Summary.R")
      # path.lst <- c("~/Desktop/HuangGroup/cvtmle_plasmode/Code/cluster/no5correct/",
      #               "~/Desktop/HuangGroup/cvtmle_plasmode/Code/cluster/no5/")
      path.lst <- c("~/Desktop/HuangGroup/cvtmle_plasmode/Code/cluster/no5/")
      situ.lst <- c("B.cor", "B")
      Effect_Size = 6.6
      mtd.lst <- c("IPW","TMLE","AIPW","DCTMLE","DCAIPW")
      factor.lst <- c("IPW","TMLE.par","TMLE.non.par","AIPW.par","AIPW.non.par",
                      "DC-TMLE.par","DC-TMLE.non.par", "DC-AIPW.par","DC-AIPW.non.par")
      violin.and.bar.plot.by.situ(path.lst, situ.lst, ci.axis.coeff=0.2)
    }
    ggsave("/Users/garethalex/Desktop/HuangGroup/cvtmle_plasmode/Output/situB.png",
           width = 14, height = 14, dpi = 300, units = "in", device='png')
  }
  
  {
    {
      source("20200816-Result-Summary.R")
      
      situ.lst <- c("C.bad")
      
      # mtd.lst <- c("GComp","TMLE","AIPW","DCTMLE","DCAIPW")
      factor.lst <- c("TMLE.par","TMLE.non.par","AIPW.par","AIPW.non.par",
                      "DC-TMLE.par","DC-TMLE.non.par", "DC-AIPW.par","DC-AIPW.non.par")
      # Effect_Size = 1.190723
      # path.lst <- c("~/Desktop/HuangGroup/cvtmle_plasmode/Code/cluster/no2/")
      # violin.and.bar.plot.by.situ(path.lst, situ.lst, ci.axis.coeff=0.2)
      Effect_Size <- 5.51788
      path.lst <- c("~/Desktop/HuangGroup/cvtmle_plasmode/Code/cluster/no2/bk/")
      # violin.and.bar.plot.by.situ(path.lst, situ.lst, ci.axis.coeff=0.2)
      violin.and.bar.plot.by.situ(path.lst, situ.lst, ci.axis.coeff=0.2,folder="ResultBefore2021-04-05")
    }
    ggsave("/Users/garethalex/Desktop/HuangGroup/cvtmle_plasmode/Output/situC.png",
           width = 14, height = 14, dpi = 300, units = "in", device='png')
  }

  plas_org <- haven::read_dta(paste0("./plas_data.dta")) 
  ncol.dat <- ncol(plas_org)
  for (i in 1:ncol.dat){
    if (is.numeric(unlist(plas_org[,i])) == F){
      print(i)
    }
  }
  covar <- plas_org[,-c(1:2,334)] # 1 is Y; 2 is A; 334 is ID
  cormat <- round(cor(covar),2)
  library(reshape2)
  melted_cormat <- melt(cormat)
  #Get the lower and upper triangles of the correlation matrix
  #Get lower triangle of the correlation matrix
  lower_tri <- cormat
  lower_tri[lower.tri(lower_tri)] <- NA #OR upper.tri function
  lower_tri
  vector.corr <- lower_tri[!is.na(lower_tri)]
  df.corr <- data.frame(vector.corr)
  ggplot(df.corr, aes(x=vector.corr)) +
    geom_histogram(color="black", fill="white")+
    scale_x_continuous(name ="Correlation")
  ggsave("/Users/garethalex/Desktop/HuangGroup/cvtmle_plasmode/Output/corr_hist.png",
         width = 10, height = 6, dpi = 300, units = "in", device='png')
  # 
  # Finished correlation matrix heatmap
  melted_cormat <- reshape2::melt(lower_tri, na.rm = TRUE)
  ggplot(data = (melted_cormat), aes(x=Var2, y=Var1, fill=(value))) + 
    geom_tile()+
    scale_fill_gradient(limit = c(-1,1), space = "Lab", 
                                     name="Pearson\nCorrelation")+
    theme(
      axis.title.x = element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.y = element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.7),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 12, barheight = 2,
                                 title.position = "top", title.hjust = 0.5))
  ggsave("/Users/garethalex/Desktop/HuangGroup/cvtmle_plasmode/Output/CovarCorr.png",
         width = 8, height = 5, dpi = 300, units = "in", device='png')
  library(tidyverse)
  thres.lst <- c(0.6,0.7,0.8,0.9)
  for (i in 1:length(thres.lst)){
    thres <- thres.lst[i]
    filt.obj <- melted_cormat %>% filter(abs(value) >thres, abs(value)!=1)
    print(nrow(filt.obj ))
  }
  filt.obj 
  
  
  ####################
  ####################
  #### Old plot codes. Not in use #####
  ####################
  ####################
  
  # {
  #   {
  #     source("20200816-Result-Summary.R")
  #     path.lst <- c("~/Desktop/HuangGroup/cvtmle_plasmode/Code/cluster/no4correct/",
  #                   "~/Desktop/HuangGroup/cvtmle_plasmode/Code/cluster/no4/",
  #                   "~/Desktop/HuangGroup/cvtmle_plasmode/Code/cluster/no4partial1/")
  #     situ.lst <- c("A.cor","A.no.int", "A.less.1st" )
  #     Effect_Size = 6.6
  #     mtd.lst <- c("IPW","TMLE","AIPW","DCTMLE","DCAIPW")
  #     violin.plot.by.situ(path.lst, situ.lst, ci.axis.coeff=0.2)
  #   }
  #   ggsave("/Users/garethalex/Desktop/HuangGroup/cvtmle_plasmode/Output/situA.png",
  #          width = 14, height = 7, dpi = 300, units = "in", device='png')
  # }
  # 
  # 
  # 
  # {
  #   {
  #     source("20200816-Result-Summary.R")
  #     path.lst <- c("~/Desktop/HuangGroup/cvtmle_plasmode/Code/cluster/no5correct/",
  #                   "~/Desktop/HuangGroup/cvtmle_plasmode/Code/cluster/no5/")
  #     situ.lst <- c("B.cor","B")
  #     Effect_Size = 6.6
  #     mtd.lst <- c("IPW","TMLE","AIPW","DCTMLE","DCAIPW")
  #     violin.plot.by.situ(path.lst, situ.lst, ci.axis.coeff=0.2)
  #   }
  #   violin.df <- violin.df.by.situ(path.lst, situ.lst, ci.axis.coeff=0.55)
  #   ggsave("/Users/garethalex/Desktop/HuangGroup/cvtmle_plasmode/Output/situB.png",
  #          width = 14, height = 7, dpi = 300, units = "in", device='png')
  # }
  # 
  # 
  # {
  #   # path.lst <- c("~/Desktop/HuangGroup/cvtmle_plasmode/Code/cluster/no2/",
  #   #               "~/Desktop/HuangGroup/cvtmle_plasmode/Code/cluster/no2partial/",
  #   #               "~/Desktop/HuangGroup/cvtmle_plasmode/Code/cluster/no2correct/")
  #   # situ.lst <- c("C.bad","C.part","C.cor")
  #   {
  #     path.lst <- c("~/Desktop/HuangGroup/cvtmle_plasmode/Code/cluster/no2/",
  #                   "~/Desktop/HuangGroup/cvtmle_plasmode/Code/cluster/no2partial/")
  #     situ.lst <- c("C.bad","C.part")
  #     Effect_Size = 1.190723
  #     mtd.lst <- c("IPW","GComp","TMLE","AIPW","DCTMLE","DCAIPW")
  #     violin.plot.by.situ(path.lst, situ.lst, ci.axis.coeff=0.55)
  #   }
  #   ggsave("/Users/garethalex/Desktop/HuangGroup/cvtmle_plasmode/Output/situC.png",
  #          width = 14, height = 7, dpi = 300, units = "in", device='png')
  # }
  # 