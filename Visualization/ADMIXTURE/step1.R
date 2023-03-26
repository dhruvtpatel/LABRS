bim=read.table("chr1-22.1KG3_phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.RAWID.filtered_2020.bim",header=F)
write.table(bim[c(2,4)], "1KGPh3_b37fwd2019_snps_b37map.txt", col.names=F,row.names=F,quote=F,sep='\t')
write.table(bim[c(2,1)], "1KGPh3_b37fwd2019_snps_b37chr.txt", col.names=F,row.names=F,quote=F,sep='\t')