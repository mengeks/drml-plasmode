##################################
## SUMMARIZE AND VISUALIZE
##################################
library(ggplot2)
library(tidyverse)
{
  setwd("~/Desktop/HuangGroup/cvtmle_plasmode/Code")
  Effect_Size <- 6.6
  source("20200816-Result-Summary.R")
  path <- "~/Desktop/HuangGroup/cvtmle_plasmode/Code/cluster/submitted092920/"
  # path <- "/Users/garethalex/Desktop/HuangGroup/cvtmle_plasmode/Code/"
  load(paste0(path,"RDataFiles/result.RData"))
  summarise.res(boot1)
}

