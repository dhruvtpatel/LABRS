library(RColorBrewer)
tbl <- read.table(​"GWAS.3.Q")
par(mar = c(1.5,4,2.5,2),cex.lab=0.75, cex.axis=0.6)

barplot(t(as.matrix(tbl)),
	col = brewer.pal(6,"Set1"), ylab = "Anc. Proportions",
	border = NA, space = 0
)
dev.off()

