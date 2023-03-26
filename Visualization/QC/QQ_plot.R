library("qqman",lib.loc="/projects/b1049/ijks/gwas/Rlibrary") 
results_log <- read.table("logistic_results.assoc_2.logistic", head=TRUE)
jpeg("QQ-Plot_logistic.jpeg")
qq(results_log$P, main = "Q-Q plot of GWAS p-values : log")
dev.off()


