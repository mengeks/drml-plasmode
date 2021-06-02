library(ggplot2)
library(tidyverse)
# setwd("~/Desktop/HuangGroup/cvtmle_plasmode/Code")
# setwd("~/Desktop/HuangGroup/cvtmle_plasmode/Code")
Effect_Size <- 6.6
# Effect_Size <- 5.51788
# Effect_Size <- 1.190723
source("~/Desktop/HuangGroup/cvtmle_plasmode/Code/20200816-Result-Summary.R")
is.timing <- T; tim.str <- ifelse(is.timing,"-timing","")
path <- "~/Desktop/HuangGroup/cvtmle_plasmode/Code/cluster/no5correct/"
(res.path <-  paste0(path,folder, "/result-","DCTMLE","-non-par",tim.str,".RData"))
load(res.path)

for (i in 1:length(boot1)){
  colnames(boot1[[i]]) <- c("ATE", "SE", "TYPE", "t", "no_cores")
}

# give.summary.res(res.path, "med_bias")
summarise.res(boot1, Effect_Size=Effect_Size)
