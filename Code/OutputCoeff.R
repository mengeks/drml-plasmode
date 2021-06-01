library(xtable)
setwd("~/Desktop/HuangGroup/cvtmle_plasmode/Code")
out.beta <- data.frame(plas_sims$TrueOutBeta[-2])
colnames(out.beta) <- "beta"
exp.beta <- data.frame(plas_sims$TrueExpBeta)
# exp.beta[12:21,] <- exp.beta[7:16,]
# exp.beta[7:11,] <- 0
exp.beta
beta.all <- cbind(out.beta, exp.beta)
colnames(beta.all) <- c("OR Coef", "PS Coef")
library(xtable)
xtable(beta.all)

library(ggplot2)
plot(out.beta[-1,1])
rownames(out.beta)
# Give var group labels and index. Include (Intercept)
# 1:"fixed", 2:"1st", 3:"int"
A.var.grp.map <- c("fixed", "1st", "int")
A.var.grp.idx <- c(rep(1, 6),rep(2, 40-4),rep(3,10))
A.var.grp <- A.var.grp.map[A.var.grp.idx]
A.var.grp <- as.factor(A.var.grp)
out.beta$A.var.grp <- A.var.grp
g0 = ggplot(out.beta[-1,], aes(x = (1:(nrow(out.beta)-1)), y = beta,color = A.var.grp))+ 
  geom_point(size = 1) + geom_hline(yintercept=0, linetype="dashed", 
                                    color = "black")+
  ggtitle("True Parameters in Scenario A, OR") +
  xlab("Variables") + ylab("Parameters") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
g0

ggsave(paste0("~/Desktop/HuangGroup/cvtmle_plasmode/Code/Plot/","A-parameter-plot.jpeg"))
