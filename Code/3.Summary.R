##################################
## SUMMARIZE AND VISUALIZE
##################################
library(ggplot2)
library(tidyverse)
{
  source("20200816-Result-Summary.R")
  # path <- "~/Desktop/HuangGroup/cvtmle_plasmode/Code/cluster/submitted091320/"
  path <- "/Users/garethalex/Desktop/HuangGroup/cvtmle_plasmode/Code/"
  load(paste0(path,"RDataFiles/091220.RData"))
  summarise.res(boot1)
}
