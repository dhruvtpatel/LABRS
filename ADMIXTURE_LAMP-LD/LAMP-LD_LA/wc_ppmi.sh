#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=5

myfiles="ppmi_gwas_preQC_chr1.vcf.gz ppmi_gwas_preQC_chr13.vcf.gz ppmi_gwas_preQC_chr17.vcf.gz ppmi_gwas_preQC_chr20.vcf.gz  ppmi_gwas_preQC_chr4.vcf.gz  ppmi_gwas_preQC_chr8.vcf.gz
ppmi_gwas_preQC_chr10.vcf.gz  ppmi_gwas_preQC_chr14.vcf.gz  ppmi_gwas_preQC_chr18.vcf.gz  ppmi_gwas_preQC_chr21.vcf.gz  ppmi_gwas_preQC_chr5.vcf.gz  ppmi_gwas_preQC_chr9.vcf.gz
ppmi_gwas_preQC_chr11.vcf.gz  ppmi_gwas_preQC_chr15.vcf.gz  ppmi_gwas_preQC_chr19.vcf.gz  ppmi_gwas_preQC_chr22.vcf.gz  ppmi_gwas_preQC_chr6.vcf.gz  ppmi_gwas_preQC_chrX.vcf.gz
ppmi_gwas_preQC_chr12.vcf.gz  ppmi_gwas_preQC_chr16.vcf.gz  ppmi_gwas_preQC_chr2.vcf.gz   ppmi_gwas_preQC_chr3.vcf.gz   ppmi_gwas_preQC_chr7.vcf.gz"

for r in $myfiles
do
zcat $r | grep -v "#" | cut -f1 | sort | uniq -c
done
